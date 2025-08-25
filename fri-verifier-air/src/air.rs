use core::borrow::Borrow;
use core::marker::PhantomData;

use p3_air::{Air, AirBuilder, BaseAir};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_matrix::Matrix;

use crate::columns::{CommitPhaseCols, num_commit_phase_cols_default};

/// AIR for FRI commit phase verification.
///
/// Verifies the FRI folding computation: f(β) = eval₀ + (β - x₀) × (eval₁ - eval₀) × (x₁ - x₀)⁻¹
#[derive(Debug)]
pub struct CommitPhaseAir<F> {
    _phantom: PhantomData<F>,
}

impl<F> CommitPhaseAir<F> {
    pub const fn new() -> Self {
        Self {
            _phantom: PhantomData,
        }
    }
}

impl<F> Default for CommitPhaseAir<F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<F: Field> BaseAir<F> for CommitPhaseAir<F> {
    fn width(&self) -> usize {
        num_commit_phase_cols_default()
    }
}

impl<AB: AirBuilder> Air<AB> for CommitPhaseAir<AB::F>
where
    AB::F: Field,
{
    #[inline]
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main
            .row_slice(0)
            .expect("main trace matrix should not be empty");
        let local = (*local).borrow();

        eval_commit_phase_constraints(builder, local);
    }
}

/// FRI folding constraint evaluation.
fn eval_commit_phase_constraints<AB: AirBuilder>(builder: &mut AB, local: &CommitPhaseCols<AB::Var, 8>)
where
    AB::Expr: PrimeCharacteristicRing,
{
    // === INDEX MANAGEMENT CONSTRAINTS ===
    // Constraint 1: parent_index = domain_index >> 1 (degree 1)
    // Note: We need to check this as a field constraint
    // For field elements, we verify: 2 * parent_index = domain_index OR 2 * parent_index + 1 = domain_index
    // This is equivalent to: domain_index - 2 * parent_index ∈ {0, 1}
    let double_parent = local.parent_index.clone() + local.parent_index.clone();
    let diff = local.domain_index.clone() - double_parent;
    let is_even = diff.clone(); // 0 if even, 1 if odd
    let is_odd = is_even.clone() - AB::Expr::ONE; // -1 if even, 0 if odd
    builder.assert_eq(is_even * is_odd, AB::Expr::ZERO); // Ensures diff ∈ {0, 1}

    // Constraint 2: sibling_index = domain_index ^ 1 (degree 1)
    // For XOR in field arithmetic: sibling = domain + (1 - 2*(domain mod 2))
    // If domain is even: sibling = domain + 1
    // If domain is odd:  sibling = domain - 1
    let sibling_expected = local.domain_index.clone() + AB::Expr::ONE - diff.clone() - diff;
    builder.assert_eq(local.sibling_index.clone(), sibling_expected);

    // === SUBGROUP POINTS CONSTRAINTS ===
    // Constraint 3: x1 = -x0 (degree 1)
    // Maps to two_adic_pcs.rs: x1 = -subgroup_start for arity-2 folding
    // For field arithmetic: x1 + x0 = 0
    builder.assert_eq(local.x1.clone() + local.x0.clone(), AB::Expr::ZERO);

    // === FRI FOLDING ARITHMETIC CONSTRAINTS ===
    // These constraints verify the actual FRI folding computation step by step
    // as implemented in two_adic_pcs.rs fold_row() method

    // Constraint 4: x_diff = x1 - x0 = -x0 - x0 = -2*x0 (degree 1)
    let expected_x_diff = -(local.x0.clone() + local.x0.clone());
    builder.assert_eq(local.x_diff.clone(), expected_x_diff);

    // Constraint 5: x_diff_inv * x_diff = 1 (degree 2)
    // Verifies that x_diff_inv is the multiplicative inverse of x_diff
    builder.assert_eq(local.x_diff_inv.clone() * local.x_diff.clone(), AB::Expr::ONE);

    // Extension field constraints (these will be verified via cross-table lookups to extension field tables)
    // For now, we document the required relationships:
    
    // Constraint 6: beta_minus_x0 = beta - x0 (extension field subtraction)
    // TODO: CTL to extension field subtraction table: (beta, x0) -> beta_minus_x0
    
    // Constraint 7: eval_diff = eval_1 - eval_0 (extension field subtraction) 
    // TODO: CTL to extension field subtraction table: (eval_1, eval_0) -> eval_diff
    
    // Constraint 8: beta_eval_product = beta_minus_x0 * eval_diff (extension field multiplication)
    // TODO: CTL to extension field multiplication table: (beta_minus_x0, eval_diff) -> beta_eval_product
    
    // Constraint 9: interpolation_term = beta_eval_product * x_diff_inv (extension field scalar multiplication)
    // TODO: CTL to extension field scalar multiplication table: (beta_eval_product, x_diff_inv) -> interpolation_term
    
    // Constraint 10: folded_eval = eval_0 + interpolation_term (extension field addition)
    // TODO: CTL to extension field addition table: (eval_0, interpolation_term) -> folded_eval
    
    // Constraint 11: beta_squared = beta * beta (extension field multiplication)
    // TODO: CTL to extension field multiplication table: (beta, beta) -> beta_squared
    
    // Constraint 12: roll_in_contribution = beta_squared * roll_in_value (extension field multiplication)
    // TODO: CTL to extension field multiplication table: (beta_squared, roll_in_value) -> roll_in_contribution

    // === DOMAIN HEIGHT CONSTRAINT ===
    // The log_height field tracks domain folding and would typically be
    // checked via cross-table lookup or public inputs
    // TODO: Add constraint linking log_height to phase progression
    
    // === INDEX BIT REPRESENTATION CONSTRAINT ===
    // The domain_index should be consistent with its bit representation used in Merkle verification
    // For a given domain_index, its binary representation should match the index_bits in MerkleTreeCols
    // This ensures consistency between FRI folding indices and Merkle proof indices
    
    // Note: The actual bit decomposition constraint will be implemented when we have
    // access to the MerkleTreeCols or a dedicated bit decomposition table
    // For now, we document the requirement that domain_index bits are properly formed

    // === FRI COMMIT PHASE MMCS VERIFICATION CONSTRAINTS ===
    // Cross-table lookup to MMCS verification table (ChallengeMmcs)
    // 
    // From fri/tests/fri.rs:791-811, the actual MMCS verification:
    // fri_mmcs.verify_batch(commit, dims, parent_index, [eval_0, sibling_eval])
    // 
    // Where:
    // - commit = fri_commit (FRI phase commitment) 
    // - dims = [width=2, height=2^log_height]
    // - parent_index = domain_index >> 1 (parent_index field)
    // - evals = [eval_0, sibling_eval] (evaluation pair)
    //
    // CTL constraint: (fri_commit, parent_index, eval_0, sibling_eval) -> mmcs_verified
    // This ensures each FRI folding step's commitment is properly verified
    
    // For now, we document the constraint structure
    // The actual CTL will be implemented when the MMCS verification table is available
    // builder.when(local.mmcs_verified.clone()).assert_eq(
    //     // CTL lookup to ChallengeMmcs table would go here
    //     AB::Expr::ONE
    // );
    
    // TODO: CTL to extension field arithmetic tables for FRI folding verification
}

/// Additional helper constraints for debugging and development.
///
/// These constraints provide sanity checks on the FRI folding arithmetic
/// and can be enabled during development to catch implementation errors.
pub fn eval_debug_constraints<AB: AirBuilder>(builder: &mut AB, local: &CommitPhaseCols<AB::Var, 8>) {
    // Verify that x1 - x0 = -2 * x0 (since x1 = -x0)
    // This is a sanity check: x1 - x0 = -x0 - x0 = -2*x0
    let expected_diff = -(local.x0.clone() + local.x0.clone());
    let actual_diff = local.x1.clone() - local.x0.clone();
    builder.assert_eq(actual_diff, expected_diff);

    // Additional debug constraints can be added here for development
}
