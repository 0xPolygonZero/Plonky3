use alloc::vec::Vec;
use p3_field::extension::BinomiallyExtendable;
use p3_field::{extension::BinomialExtensionField, BasedVectorSpace, Field, PrimeCharacteristicRing};
use p3_matrix::dense::RowMajorMatrix;

use crate::columns::CommitPhaseCols;

/// Data structure representing a single FRI commit phase folding step,
/// aligned with `CommitPhaseCols`.
#[derive(Debug, Clone)]
pub struct CommitPhaseStep<F: Field> {
    // Identification (base)
    pub query_index: F,
    pub phase_index: F,

    // Challenges & evals (ext)
    pub beta: BinomialExtensionField<F, 4>,
    pub eval_0: BinomialExtensionField<F, 4>,
    pub eval_1: BinomialExtensionField<F, 4>,

    // Roll-in gating & values
    pub is_roll_in: F,                                 // {0,1}
    pub reduced_opening: BinomialExtensionField<F, 4>, // ro at this height (or 0)
    pub beta_squared: BinomialExtensionField<F, 4>,
    pub roll_in_contribution: BinomialExtensionField<F, 4>,

    // Indices (base)
    pub domain_index: F,
    pub parent_index: F,
    pub sibling_index: F,
    pub lsb: F, // domain_index & 1

    // Height (base)
    pub log_height: F,

    // Subgroup points (base)
    pub x0: F,
    pub x1: F,

    // Folding intermediates
    pub beta_minus_x0: BinomialExtensionField<F, 4>,
    pub eval_diff: BinomialExtensionField<F, 4>,
    pub x_diff: F,
    pub x_diff_inv: F,
    pub beta_eval_product: BinomialExtensionField<F, 4>,
    pub interpolation_term: BinomialExtensionField<F, 4>,
    pub folded_eval: BinomialExtensionField<F, 4>,

    // Commitment digest for this phase
    pub fri_commit: [F; 8],
}

/// Generate trace rows for the CommitPhase table (aligned with `CommitPhaseCols`).
pub fn generate_commit_phase_trace<F>(steps: &[CommitPhaseStep<F>]) -> RowMajorMatrix<F>
where
    F: Field + BinomiallyExtendable<4>,
{
    let num_rows = steps.len();
    let num_cols = core::mem::size_of::<CommitPhaseCols<u8, 8>>();

    let mut trace_data = Vec::with_capacity(num_rows * num_cols);

    for step in steps {
        // Ext components
        let beta_c = step.beta.as_basis_coefficients_slice();
        let eval0_c = step.eval_0.as_basis_coefficients_slice();
        let eval1_c = step.eval_1.as_basis_coefficients_slice();
        let ro_c = step.reduced_opening.as_basis_coefficients_slice();
        let beta2_c = step.beta_squared.as_basis_coefficients_slice();
        let rollin_c = step.roll_in_contribution.as_basis_coefficients_slice();
        let bmx0_c = step.beta_minus_x0.as_basis_coefficients_slice();
        let ediff_c = step.eval_diff.as_basis_coefficients_slice();
        let beprod_c = step.beta_eval_product.as_basis_coefficients_slice();
        let interp_c = step.interpolation_term.as_basis_coefficients_slice();
        let folded_c = step.folded_eval.as_basis_coefficients_slice();

        // === Order must match CommitPhaseCols ===
        // Identifiers
        trace_data.push(step.query_index);
        trace_data.push(step.phase_index);

        // beta, eval_0, eval_1
        trace_data.extend_from_slice(beta_c);
        trace_data.extend_from_slice(eval0_c);
        trace_data.extend_from_slice(eval1_c);

        // roll-in gate & values
        trace_data.push(step.is_roll_in);
        trace_data.extend_from_slice(ro_c);
        trace_data.extend_from_slice(beta2_c);
        trace_data.extend_from_slice(rollin_c);

        // indices
        trace_data.push(step.domain_index);
        trace_data.push(step.parent_index);
        trace_data.push(step.sibling_index);
        trace_data.push(step.lsb);

        // height
        trace_data.push(step.log_height);

        // x-points
        trace_data.push(step.x0);
        trace_data.push(step.x1);

        // intermediates
        trace_data.extend_from_slice(bmx0_c);
        trace_data.extend_from_slice(ediff_c);
        trace_data.push(step.x_diff);
        trace_data.push(step.x_diff_inv);
        trace_data.extend_from_slice(beprod_c);
        trace_data.extend_from_slice(interp_c);
        trace_data.extend_from_slice(folded_c);

        // commitment digest
        trace_data.extend_from_slice(&step.fri_commit);
    }

    RowMajorMatrix::new(trace_data, num_cols)
}

/// Parameters for constructing a `CommitPhaseStep`.
pub struct CommitPhaseStepParams<F: Field> {
    pub query_index: usize,
    pub phase_index: usize,

    pub beta: BinomialExtensionField<F, 4>,
    pub eval_0: BinomialExtensionField<F, 4>,
    pub eval_1: BinomialExtensionField<F, 4>,

    pub reduced_opening: BinomialExtensionField<F, 4>, // ro at this height (0 if none)
    pub is_roll_in: bool,

    pub domain_index: usize,
    pub log_height: usize,
    pub x0: F, // subgroup start for the row

    // Intermediates (you may precompute externally if desired)
    pub beta_minus_x0: BinomialExtensionField<F, 4>,
    pub eval_diff: BinomialExtensionField<F, 4>,
    pub x_diff: F,
    pub x_diff_inv: F,
    pub beta_eval_product: BinomialExtensionField<F, 4>,
    pub interpolation_term: BinomialExtensionField<F, 4>,
    pub folded_eval_pre_rollin: BinomialExtensionField<F, 4>, // before adding roll-in

    pub fri_commit: [F; 8],
}

/// Construct a `CommitPhaseStep` and compute β² and roll-in contribution.
pub fn create_commit_phase_step<F>(p: CommitPhaseStepParams<F>) -> CommitPhaseStep<F>
where
    F: Field + BinomiallyExtendable<4>,
{
    let beta_squared = p.beta * p.beta;
    let roll_in_contribution = if p.is_roll_in {
        beta_squared * p.reduced_opening
    } else {
        BinomialExtensionField::ZERO
    };

    let folded_eval = p.folded_eval_pre_rollin + roll_in_contribution;

    let di = p.domain_index;
    let parent = di >> 1;
    let lsb_u = di & 1;

    CommitPhaseStep {
        query_index: F::from_usize(p.query_index),
        phase_index: F::from_usize(p.phase_index),

        beta: p.beta,
        eval_0: p.eval_0,
        eval_1: p.eval_1,

        is_roll_in: if p.is_roll_in { F::ONE } else { F::ZERO },
        reduced_opening: p.reduced_opening,
        beta_squared,
        roll_in_contribution,

        domain_index: F::from_usize(di),
        parent_index: F::from_usize(parent),
        sibling_index: F::from_usize(di ^ 1),
        lsb: F::from_usize(lsb_u),

        log_height: F::from_usize(p.log_height),

        x0: p.x0,
        x1: -p.x0,

        beta_minus_x0: p.beta_minus_x0,
        eval_diff: p.eval_diff,
        x_diff: p.x_diff,
        x_diff_inv: p.x_diff_inv,
        beta_eval_product: p.beta_eval_product,
        interpolation_term: p.interpolation_term,
        folded_eval,

        fri_commit: p.fri_commit,
    }
}

/// Extract `CommitPhaseStep`s from a real FRI proof, binding to the transcript
/// in the same order as verification. Assumes `reduced_openings` is provided as
/// `(log_height, ro)` pairs sorted **descending** by height starting at `log_max_height`.
pub fn extract_commit_phase_steps_from_fri_proof<F, M, C, Witness, InputProof>(
    fri_proof: &p3_fri::FriProof<F, M, Witness, InputProof>,
    challenger: &mut C,
    fri_params: &p3_fri::FriParameters<M>,
    query_slot: usize,
    reduced_openings: &[(usize, BinomialExtensionField<F, 4>)],
) -> Vec<CommitPhaseStep<F>>
where
    F: Field
        + BinomiallyExtendable<4>
        + p3_field::TwoAdicField
        + p3_field::PrimeCharacteristicRing,
    M: p3_commit::Mmcs<F>,
    M::Commitment: AsRef<[F; 8]>,
    C: p3_challenger::FieldChallenger<F>
        + p3_challenger::CanObserve<M::Commitment>
        + p3_challenger::CanSampleBits<usize>,
    Witness: Clone,
{
    let mut steps = Vec::new();

    // Bind betas: observe commits then sample βᵢ
    let betas: Vec<BinomialExtensionField<F, 4>> = fri_proof
        .commit_phase_commits
        .iter()
        .map(|c| {
            challenger.observe(c.clone());
            challenger.sample_algebra_element()
        })
        .collect();

    // Observe final polynomial coefficients (binding the transcript)
    for coeff in &fri_proof.final_poly {
        challenger.observe_algebra_element(*coeff);
    }

    // Sample the random index (log_max_height bits) exactly like the verifier
    let log_max_height =
        fri_proof.commit_phase_commits.len() + fri_params.log_blowup + fri_params.log_final_poly_len;
    let mut domain_index = challenger.sample_bits(log_max_height);

    // Locate the query proof we’re extracting (e.g., 0 for single-query setups)
    let Some(qp) = fri_proof.query_proofs.get(query_slot) else {
        return steps;
    };

    // Initialize folded_eval from the reduced opening at max height
    let mut ro_iter = reduced_openings.iter().copied().peekable();
    let mut current_folded =
        if let Some((lh, ro)) = ro_iter.next() {
            debug_assert_eq!(lh, log_max_height, "First reduced opening must be at max height");
            ro
        } else {
            // No reduced opening provided; treat as zero (shouldn’t happen in practice).
            BinomialExtensionField::ZERO
        };

    for (phase_idx, (commit, opening)) in fri_proof
        .commit_phase_commits
        .iter()
        .zip(qp.commit_phase_openings.iter())
        .enumerate()
    {
        let beta = betas[phase_idx];

        // Height and subgroup points for this phase (arity = 2)
        let log_folded_height = log_max_height - phase_idx - 1;
        let rev_bits = p3_util::reverse_bits_len(domain_index, log_folded_height + 1);
        let generator = F::two_adic_generator(log_folded_height + 1);
        let x0 = generator.exp_u64(rev_bits as u64);
        let x1 = -x0;

        // eval_0 is the current folded_eval; eval_1 is sibling from proof
        let eval_0 = current_folded;
        let eval_1: BinomialExtensionField<F, 4> = opening.sibling_value.into();

        // β - x0
        let x0_ext =
            BinomialExtensionField::from_basis_coefficients_slice(&[x0, F::ZERO, F::ZERO, F::ZERO])
                .expect("valid x0 ext");
        let beta_minus_x0 = beta - x0_ext;

        // eval_1 - eval_0
        let eval_diff = eval_1 - eval_0;

        // x1 - x0 and its inverse
        let x_diff = x1 - x0;
        let x_diff_inv = x_diff.inverse();

        // (β - x0) * (eval_1 - eval_0)
        let beta_eval_product = beta_minus_x0 * eval_diff;

        // interpolation_term = ((β - x0)*(eval_1 - eval_0)) * (x1 - x0)^(-1)
        let x_diff_inv_ext = BinomialExtensionField::from_basis_coefficients_slice(&[
            x_diff_inv,
            F::ZERO,
            F::ZERO,
            F::ZERO,
        ])
        .expect("valid inv ext");
        let interpolation_term = beta_eval_product * x_diff_inv_ext;

        // folded_eval before roll-in
        let folded_eval_pre_rollin = eval_0 + interpolation_term;

        // β² and roll-in at this height (if any)
        let beta_squared = beta * beta;
        let (is_roll_in, reduced_opening, roll_in_contribution) =
            if let Some((_lh, ro)) = ro_iter.next_if(|(lh, _)| *lh == log_folded_height) {
                (F::ONE, ro, beta_squared * ro)
            } else {
                (F::ZERO, BinomialExtensionField::ZERO, BinomialExtensionField::ZERO)
            };

        // Update folded value for the next phase
        current_folded = folded_eval_pre_rollin + roll_in_contribution;

        // Commit digest
        let fri_commit: [F; 8] = *commit.as_ref();

        // Build step
        let step = CommitPhaseStep {
            query_index: F::from_usize(query_slot),
            phase_index: F::from_usize(phase_idx),

            beta,
            eval_0,
            eval_1,

            is_roll_in,
            reduced_opening,
            beta_squared,
            roll_in_contribution,

            domain_index: F::from_usize(domain_index),
            parent_index: F::from_usize(domain_index >> 1),
            sibling_index: F::from_usize(domain_index ^ 1),
            lsb: F::from_usize(domain_index & 1),

            log_height: F::from_usize(log_folded_height),

            x0,
            x1,

            beta_minus_x0,
            eval_diff,
            x_diff,
            x_diff_inv,
            beta_eval_product,
            interpolation_term,
            folded_eval: current_folded,

            fri_commit,
        };
        steps.push(step);

        // Move to parent index for next round
        domain_index >>= 1;
    }

    steps
}
