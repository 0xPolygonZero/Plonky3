use core::borrow::Borrow;
use core::marker::PhantomData;

use p3_air::{Air, AirBuilder, BaseAir};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_matrix::Matrix;

use crate::columns::{CommitPhaseCols, num_commit_phase_cols};

/// AIR for the CommitPhase table of FRI verifier.
///
/// Each row represents one FRI folding step in the commit phase.
/// This implements the constraints that verify the FRI folding arithmetic
/// as performed in simulate_fri_verification_immediate_values.
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

impl<F: Field> BaseAir<F> for CommitPhaseAir<F> {
    fn width(&self) -> usize {
        num_commit_phase_cols()
    }
}

impl<AB: AirBuilder> Air<AB> for CommitPhaseAir<AB::F>
where
    AB::F: Field,
{
    #[inline]
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main.row_slice(0).expect("The matrix is empty?");
        let local = (*local).borrow();

        eval_commit_phase_constraints(builder, local);
    }
}

/// Evaluates the constraints for the CommitPhase AIR.
///
/// Maps directly to the FRI folding logic from fri.rs:727-908
fn eval_commit_phase_constraints<AB: AirBuilder>(builder: &mut AB, local: &CommitPhaseCols<AB::Var>)
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
    // Maps to fri.rs:848: x1 = -subgroup_start
    // For field arithmetic: x1 + x0 = 0
    builder.assert_eq(local.x1.clone() + local.x0.clone(), AB::Expr::ZERO);

    // === INTERMEDIATE COMPUTATION CONSTRAINTS ===
    // These break down the FRI folding formula into degree ≤ 3 constraints
    // Maps to fri.rs:857-874: folded_eval = e0 + (beta - x0) * (e1 - e0) * (x1 - x0)^(-1)

    // Constraint 4: beta_minus_x0 = beta - x0 (degree 1)
    builder.assert_eq(
        local.beta_minus_x0.clone(),
        local.beta.clone() - local.x0.clone(),
    );

    // Constraint 5: eval_1_minus_eval_0 = eval_1 - eval_0 (degree 1)
    builder.assert_eq(
        local.eval_1_minus_eval_0.clone(),
        local.eval_1.clone() - local.eval_0.clone(),
    );

    // Constraint 6: x1_minus_x0 = x1 - x0 (degree 1)
    builder.assert_eq(
        local.x1_minus_x0.clone(),
        local.x1.clone() - local.x0.clone(),
    );

    // Constraint 7: inverse * x1_minus_x0 = 1 (degree 2)
    // This verifies that inverse = (x1 - x0)^(-1)
    builder.assert_eq(
        local.inverse.clone() * local.x1_minus_x0.clone(),
        AB::Expr::ONE,
    );

    // Constraint 8: intermediate = beta_minus_x0 * eval_1_minus_eval_0 * inverse (degree 3)
    // This is the core FRI folding computation
    builder.assert_eq(
        local.intermediate.clone(),
        local.beta_minus_x0.clone().into()
            * local.eval_1_minus_eval_0.clone().into()
            * local.inverse.clone().into(),
    );

    // Constraint 9: folded_result_pre = eval_0 + intermediate (degree 1)
    builder.assert_eq(
        local.folded_result_pre.clone(),
        local.eval_0.clone() + local.intermediate.clone(),
    );

    // === ROLL-IN CONSTRAINT ===
    // Constraint 10: folded_result = folded_result_pre + roll_in_value (degree 1)
    // Maps to fri.rs:898: folded_eval += roll_in
    builder.assert_eq(
        local.folded_result.clone(),
        local.folded_result_pre.clone() + local.roll_in_value.clone(),
    );

    // === DOMAIN HEIGHT CONSTRAINT ===
    // Constraint 11: Verify log_height is consistent
    // This would typically be checked via cross-table lookup or public inputs
    // For now, we can add range constraints or consistency checks

    // Note: MMCS verification constraint is intentionally omitted as it would
    // require complex cryptographic constraints. In practice, this would be
    // verified via witness generation and potentially cross-table lookups.
}

/// Additional helper constraints that might be useful for debugging
pub fn eval_debug_constraints<AB: AirBuilder>(builder: &mut AB, local: &CommitPhaseCols<AB::Var>) {
    // Verify that x1_minus_x0 = -2 * x0 (since x1 = -x0)
    // This is a sanity check: x1 - x0 = -x0 - x0 = -2*x0
    let expected_diff = -(local.x0.clone() + local.x0.clone());
    builder.assert_eq(local.x1_minus_x0.clone(), expected_diff);

    // Verify that sibling_index is actually the XOR of domain_index with 1
    // This is a more explicit check of the XOR constraint
}
