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

impl<F> Default for CommitPhaseAir<F> {
    fn default() -> Self {
        Self::new()
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
        let local = main
            .row_slice(0)
            .expect("main trace matrix should not be empty");
        let local = (*local).borrow();

        eval_commit_phase_constraints(builder, local);
    }
}

/// Evaluates the constraints for the CommitPhase AIR.
///
/// Maps directly to the FRI folding logic from verifier.rs and two_adic_pcs.rs
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
    // Maps to two_adic_pcs.rs: x1 = -subgroup_start for arity-2 folding
    // For field arithmetic: x1 + x0 = 0
    builder.assert_eq(local.x1.clone() + local.x0.clone(), AB::Expr::ZERO);

    // === EXTENSION FIELD ELEMENT CONSTRAINTS ===
    // Extension field arithmetic operations (beta - x0, eval_1 - eval_0, etc.)
    // will be handled by dedicated extension field operation tables.
    // This table focuses on the high-level FRI folding logic and uses
    // cross-table lookups to verify the extension field computations.

    // For now, we don't have explicit FRI folding constraints here since
    // the extension field operations will be verified via cross-table lookups
    // to dedicated extension field arithmetic tables.

    // === DOMAIN HEIGHT CONSTRAINT ===
    // The log_height field tracks domain folding and would typically be
    // checked via cross-table lookup or public inputs
    // TODO: Add constraint linking log_height to phase progression

    // TODO: CTL to MMCS verification table
    // TODO: CTL to extension field arithmetic tables for FRI folding verification
}

/// Additional helper constraints that might be useful for debugging
pub fn eval_debug_constraints<AB: AirBuilder>(builder: &mut AB, local: &CommitPhaseCols<AB::Var>) {
    // Verify that x1 - x0 = -2 * x0 (since x1 = -x0)
    // This is a sanity check: x1 - x0 = -x0 - x0 = -2*x0
    let expected_diff = -(local.x0.clone() + local.x0.clone());
    let actual_diff = local.x1.clone() - local.x0.clone();
    builder.assert_eq(actual_diff, expected_diff);

    // Additional debug constraints can be added here for development
}
