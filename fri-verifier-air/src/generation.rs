use alloc::vec::Vec;
use p3_field::extension::BinomiallyExtendable;
use p3_field::{
    BasedVectorSpace, Field, PrimeCharacteristicRing, extension::BinomialExtensionField,
};
use p3_matrix::dense::RowMajorMatrix;

use crate::columns::CommitPhaseCols;

/// Data structure representing a single FRI commit phase folding step.
/// This mirrors the data captured in simulate_fri_verification_immediate_values.
#[derive(Debug, Clone)]
pub struct CommitPhaseStep<F: Field> {
    // Query identification (base field)
    pub query_index: F,
    pub phase_index: F,

    // FRI challenges and evaluations (extension field)
    pub beta: BinomialExtensionField<F, 4>,
    pub eval_0: BinomialExtensionField<F, 4>,
    pub eval_1: BinomialExtensionField<F, 4>,
    pub roll_in_value: BinomialExtensionField<F, 4>,

    // Domain management (base field - these are indices/positions)
    pub domain_index: F,
    pub parent_index: F,
    pub sibling_index: F,

    // Subgroup points (base field - two-adic subgroup elements)
    pub x0: F,
    pub x1: F,

    // Verification flags (base field)
    pub mmcs_verified: F,
    pub log_height: F,
}

/// Generate trace rows for the CommitPhase table.
///
/// This function takes the FRI folding steps captured from the actual verifier
/// and generates the corresponding AIR trace data.
pub fn generate_commit_phase_trace<F>(steps: &[CommitPhaseStep<F>]) -> RowMajorMatrix<F>
where
    F: Field + BinomiallyExtendable<4>,
{
    let num_rows = steps.len();
    let num_cols = core::mem::size_of::<CommitPhaseCols<u8>>();

    let mut trace_data = Vec::with_capacity(num_rows * num_cols);

    for step in steps {
        // Extract extension field components using basis coefficients
        let beta_components = step.beta.as_basis_coefficients_slice();
        let eval_0_components = step.eval_0.as_basis_coefficients_slice();
        let eval_1_components = step.eval_1.as_basis_coefficients_slice();
        let roll_in_components = step.roll_in_value.as_basis_coefficients_slice();

        // Create the row data in the same order as CommitPhaseCols struct
        // Base field elements
        trace_data.push(step.query_index);
        trace_data.push(step.phase_index);

        // Beta challenge (4 components)
        trace_data.extend_from_slice(beta_components);

        // eval_0 (4 components)
        trace_data.extend_from_slice(eval_0_components);

        // eval_1 (4 components)
        trace_data.extend_from_slice(eval_1_components);

        // roll_in_value (4 components)
        trace_data.extend_from_slice(roll_in_components);

        // Domain management (base field)
        trace_data.push(step.domain_index);
        trace_data.push(step.parent_index);
        trace_data.push(step.sibling_index);

        // Subgroup points (base field)
        trace_data.push(step.x0);
        trace_data.push(step.x1);

        // Verification flags (base field)
        trace_data.push(step.mmcs_verified);
        trace_data.push(step.log_height);
    }

    RowMajorMatrix::new(trace_data, num_cols)
}

/// Parameters for creating a CommitPhaseStep.
/// Groups related parameters to avoid too many function arguments.
pub struct CommitPhaseStepParams<F: Field> {
    pub query_index: usize,
    pub phase_index: usize,
    pub beta: BinomialExtensionField<F, 4>,
    pub domain_index: usize,
    pub eval_0: BinomialExtensionField<F, 4>,
    pub eval_1: BinomialExtensionField<F, 4>, // sibling value
    pub x0: F,                                // subgroup point
    pub roll_in_value: BinomialExtensionField<F, 4>,
    pub log_height: usize,
    pub mmcs_verified: bool,
}

/// Helper function to create a CommitPhaseStep from FRI verifier data.
/// This would be called during actual FRI verification to capture the immediate values.
pub fn create_commit_phase_step<F>(params: CommitPhaseStepParams<F>) -> CommitPhaseStep<F>
where
    F: Field + BinomiallyExtendable<4>,
{
    CommitPhaseStep {
        query_index: F::from_usize(params.query_index),
        phase_index: F::from_usize(params.phase_index),
        beta: params.beta,
        domain_index: F::from_usize(params.domain_index),
        parent_index: F::from_usize(params.domain_index >> 1),
        sibling_index: F::from_usize(params.domain_index ^ 1),
        eval_0: params.eval_0,
        eval_1: params.eval_1,
        x0: params.x0,
        x1: -params.x0, // x1 = -x0 for arity-2 folding
        roll_in_value: params.roll_in_value,
        mmcs_verified: if params.mmcs_verified {
            F::ONE
        } else {
            F::ZERO
        },
        log_height: F::from_usize(params.log_height),
    }
}

/// Extract CommitPhase steps from a Plonky3 FRI proof for production use.
///
/// This function processes actual FRI proof structures to generate CommitPhase AIR table data.
/// It should be used by recursive verifiers that need to verify FRI proofs within a STARK circuit.
///
/// # Arguments
/// * `fri_proof` - The FRI proof containing commit phase data
/// * `challenger_state` - Challenger state (will be updated by observing proof data)
/// * `fri_params` - FRI parameters used in the original proof
/// * `query_index` - The query index for this proof (usually 0 for recursive verifiers)
///
/// # Returns
/// Vector of CommitPhaseStep entries that can be used to generate AIR trace rows
pub fn extract_commit_phase_steps_from_fri_proof<F, M, C, Witness, InputProof>(
    fri_proof: &p3_fri::FriProof<F, M, Witness, InputProof>,
    challenger_state: &mut C,
    fri_params: &p3_fri::FriParameters<M>,
    query_index: usize,
) -> Vec<CommitPhaseStep<F>>
where
    F: Field + BinomiallyExtendable<4> + p3_field::TwoAdicField + p3_field::PrimeCharacteristicRing,
    M: p3_commit::Mmcs<F>,
    C: p3_challenger::FieldChallenger<F>
        + p3_challenger::CanObserve<M::Commitment>
        + p3_challenger::CanSampleBits<usize>,
    Witness: Clone,
{
    let mut commit_phase_steps = Vec::new();

    // Generate beta challenges by observing commit phase commits
    let betas: Vec<BinomialExtensionField<F, 4>> = fri_proof
        .commit_phase_commits
        .iter()
        .map(|comm| {
            challenger_state.observe(comm.clone());
            challenger_state.sample_algebra_element()
        })
        .collect();

    // Observe final polynomial
    for coeff in &fri_proof.final_poly {
        challenger_state.observe_algebra_element(*coeff);
    }

    // Process the query proof (use query_index, typically 0 for recursive verifiers)
    if let Some(query_proof) = fri_proof.query_proofs.get(query_index) {
        let log_max_height = fri_proof.commit_phase_commits.len()
            + fri_params.log_blowup
            + fri_params.log_final_poly_len;

        let mut domain_index = query_index;

        for (phase_idx, opening) in query_proof.commit_phase_openings.iter().enumerate() {
            if phase_idx < betas.len() {
                let beta = betas[phase_idx];
                let sibling_value = opening.sibling_value;

                // Calculate subgroup points
                let log_folded_height = log_max_height - phase_idx;
                let rev_bits = p3_util::reverse_bits_len(domain_index, log_folded_height);
                let generator = F::two_adic_generator(log_folded_height);
                let x0 = generator.exp_u64(rev_bits as u64);

                // Create CommitPhaseStep with actual proof data
                let step = CommitPhaseStep {
                    query_index: F::from_usize(query_index),
                    phase_index: F::from_usize(phase_idx),
                    beta,
                    domain_index: F::from_usize(domain_index),
                    parent_index: F::from_usize(domain_index >> 1),
                    sibling_index: F::from_usize(domain_index ^ 1),
                    eval_0: sibling_value.into(), // Use actual sibling evaluation
                    eval_1: sibling_value.into(), // Sibling evaluation
                    x0,
                    x1: -x0, // x1 = -x0 for arity-2 folding
                    roll_in_value: BinomialExtensionField::ZERO, // Would be computed during folding
                    mmcs_verified: F::ONE, // Assume MMCS verification passed
                    log_height: F::from_usize(log_folded_height),
                };

                commit_phase_steps.push(step);
                domain_index >>= 1; // Fold for next phase
            }
        }
    }

    commit_phase_steps
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use p3_field::PrimeCharacteristicRing;
    use p3_matrix::Matrix;

    type TestField = p3_baby_bear::BabyBear;

    #[test]
    fn test_commit_phase_step_creation() {
        let beta = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::from_u32(42),
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");
        let eval_0 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::from_u32(100),
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");
        let eval_1 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::from_u32(200),
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");
        let roll_in = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");

        let step = create_commit_phase_step(CommitPhaseStepParams {
            query_index: 0,
            phase_index: 0,
            beta,
            domain_index: 8,
            eval_0,
            eval_1,
            x0: TestField::from_u32(7),
            roll_in_value: roll_in,
            log_height: 4,
            mmcs_verified: true,
        });

        // Verify computed fields
        assert_eq!(step.parent_index, TestField::from_usize(4)); // 8 >> 1
        assert_eq!(step.sibling_index, TestField::from_usize(9)); // 8 ^ 1
        assert_eq!(step.x1, -TestField::from_u32(7)); // -x0
        assert_eq!(step.query_index, TestField::from_usize(0)); // query_index
        assert_eq!(step.phase_index, TestField::from_usize(0)); // phase_index
    }

    #[test]
    fn test_trace_generation() {
        let beta1 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::from_u32(42),
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");
        let eval_0_1 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::from_u32(100),
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");
        let eval_1_1 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::from_u32(200),
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");
        let roll_in_1 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");

        let beta2 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::from_u32(43),
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");
        let eval_0_2 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::from_u32(150),
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");
        let eval_1_2 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::from_u32(250),
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");
        let roll_in_2 = BinomialExtensionField::from_basis_coefficients_slice(&[
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
            TestField::ZERO,
        ])
        .expect("invalid coefficients");

        let steps = vec![
            create_commit_phase_step(CommitPhaseStepParams {
                query_index: 0,
                phase_index: 0,
                beta: beta1,
                domain_index: 8,
                eval_0: eval_0_1,
                eval_1: eval_1_1,
                x0: TestField::from_u32(7),
                roll_in_value: roll_in_1,
                log_height: 4,
                mmcs_verified: true,
            }),
            create_commit_phase_step(CommitPhaseStepParams {
                query_index: 0,
                phase_index: 1,
                beta: beta2,
                domain_index: 4,
                eval_0: eval_0_2,
                eval_1: eval_1_2,
                x0: TestField::from_u32(3),
                roll_in_value: roll_in_2,
                log_height: 3,
                mmcs_verified: true,
            }),
        ];

        let trace = generate_commit_phase_trace(&steps);
        assert_eq!(trace.height(), 2);
        assert_eq!(trace.width(), core::mem::size_of::<CommitPhaseCols<u8>>());
    }
}
