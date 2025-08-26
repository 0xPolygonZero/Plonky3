use core::borrow::Borrow;
use core::marker::PhantomData;

use p3_air::{Air, AirBuilder, BaseAir};
use p3_field::{Field, PrimeCharacteristicRing};
use p3_matrix::Matrix;

use crate::columns::{num_commit_phase_cols_default, CommitPhaseCols};

/// AIR for FRI commit phase verification.
///
/// Verifies one arity-2 FRI folding step per row:
///   folded_eval = eval₀
///                + ((β - x₀) * (eval₁ - eval₀)) * (x₁ - x₀)⁻¹
///                + is_roll_in * (β² * reduced_opening)
///
/// Notes:
/// - Index relations are enforced with an explicit `lsb` bit:
///     domain_index = 2 * parent_index + lsb
///     sibling_index = domain_index + (1 - 2*lsb)
/// - `eval_1` is the sibling evaluation proven against `fri_commit` via an MMCS chip (CTL).
/// - All extension-field ops over 4-limb representations are expected to be enforced
///   via CTLs to extension arithmetic chips (add/sub/mul/scale).
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
        // (Optional) extra sanity checks during development:
        // eval_debug_constraints(builder, local);
    }
}

/// FRI folding constraint evaluation for one row.
fn eval_commit_phase_constraints<AB: AirBuilder>(
    builder: &mut AB,
    local: &CommitPhaseCols<AB::Var, 8>,
) where
    AB::Expr: PrimeCharacteristicRing,
{
    // =========================
    // INDEX MANAGEMENT
    // =========================
    // lsb ∈ {0,1}
    builder.assert_eq(local.lsb.clone() * (local.lsb.clone() - AB::Expr::ONE), AB::Expr::ZERO);

    // domain_index = 2 * parent_index + lsb
    let two_parent = local.parent_index.clone() + local.parent_index.clone();
    builder.assert_eq(local.domain_index.clone(), two_parent.clone() + local.lsb.clone().into());

    // sibling_index = domain_index + (1 - 2*lsb)
    let two_lsb = local.lsb.clone() + local.lsb.clone();
    let sibling_expected = local.domain_index.clone() + AB::Expr::ONE - two_lsb;
    builder.assert_eq(local.sibling_index.clone(), sibling_expected);

    // =========================
    // SUBGROUP POINTS
    // =========================
    // x1 = -x0
    builder.assert_eq(local.x1.clone() + local.x0.clone(), AB::Expr::ZERO);

    // x_diff = x1 - x0 = -2 * x0
    let expected_x_diff = -(local.x0.clone() + local.x0.clone());
    builder.assert_eq(local.x_diff.clone(), expected_x_diff);

    // x_diff_inv * x_diff = 1
    builder.assert_eq(local.x_diff_inv.clone() * local.x_diff.clone(), AB::Expr::ONE);

    // (Optionally enforce two-adic subgroup membership for x0 via a dedicated chip/CTL):
    // TODO(CTL): assert x0^(2^(log_height+1)) = 1 and x0 != 0

    // =========================
    // EXTENSION-FIELD ARITHMETIC (4-limb rep)
    // All below should be wired with CTLs into extension-arith chips.
    // We list the identities that must hold.
    // =========================

    // beta_minus_x0 = beta - x0
    // TODO(CTL): (beta, x0) -> beta_minus_x0

    // eval_diff = eval_1 - eval_0
    // TODO(CTL): (eval_1, eval_0) -> eval_diff

    // beta_eval_product = beta_minus_x0 * eval_diff
    // TODO(CTL): (beta_minus_x0, eval_diff) -> beta_eval_product

    // interpolation_term = beta_eval_product * x_diff_inv
    // (scalar mult by base field element)
    // TODO(CTL): (beta_eval_product, x_diff_inv) -> interpolation_term

    // beta_squared = beta * beta
    // TODO(CTL): (beta, beta) -> beta_squared

    // roll_in_contribution = is_roll_in * (beta_squared * reduced_opening)
    // Break into two steps to enable CTL reuse:
    //   tmp = beta_squared * reduced_opening
    //   roll_in_contribution = is_roll_in * tmp
    // TODO(CTL): (beta_squared, reduced_opening) -> tmp
    // TODO(CTL): gate_by_bool(is_roll_in, tmp) -> roll_in_contribution

    // folded_eval = eval_0 + interpolation_term + roll_in_contribution
    // We can enforce as two chained additions:
    //   tmp2 = eval_0 + interpolation_term
    //   folded_eval = tmp2 + roll_in_contribution
    // TODO(CTL): (eval_0, interpolation_term) -> tmp2
    // TODO(CTL): (tmp2, roll_in_contribution) -> folded_eval

    // =========================
    // HEIGHT / PHASE PROGRESSION
    // =========================
    // log_height should follow: log_height = log_max_height - phase_index - 1
    // Usually enforced via public inputs / per-phase meta or a dedicated chip.
    // TODO(CTL): tie (phase_index) to (log_height) per query_index.

    // =========================
    // MMCS VERIFICATION (FRI COMMIT PHASE)
    // =========================
    // Each row must be backed by an MMCS proof:
    //   fri_mmcs.verify_batch(
    //       commit = fri_commit,
    //       dims   = [width=2, height=2^log_height],
    //       index  = parent_index,
    //       evals  = [eval_0, eval_1],
    //       proof  = opening_proof (lives in MMCS chip)
    //   )
    //
    // Enforce via a CTL to the MMCS chip keyed by (query_index, phase_index):
    // TODO(CTL): (fri_commit, log_height, parent_index, [eval_0, eval_1]) ∈ MMCS.verify
}

/// Optional sanity checks for development.
#[allow(dead_code)]
pub fn eval_debug_constraints<AB: AirBuilder>(
    builder: &mut AB,
    local: &CommitPhaseCols<AB::Var, 8>,
) {
    // x1 - x0 = -2*x0
    let expected_diff = -(local.x0.clone() + local.x0.clone());
    let actual_diff = local.x1.clone() - local.x0.clone();
    builder.assert_eq(actual_diff, expected_diff);

    // lsb boolean
    builder.assert_eq(local.lsb.clone() * (local.lsb.clone() - AB::Expr::ONE), AB::Expr::ZERO);

    // domain_index relation
    let two_parent = local.parent_index.clone() + local.parent_index.clone();
    builder.assert_eq(local.domain_index.clone(), two_parent + local.lsb.clone());
}
