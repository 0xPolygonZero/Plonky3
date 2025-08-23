use alloc::vec::Vec;
use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;

use crate::columns::CommitPhaseCols;

/// Data structure representing a single FRI commit phase folding step.
/// This mirrors the data captured in simulate_fri_verification_immediate_values.
#[derive(Debug, Clone)]
pub struct CommitPhaseStep<F: Field> {
    pub query_id: F,
    pub phase_idx: F,
    pub beta: F,
    pub domain_index: F,
    pub parent_index: F,
    pub sibling_index: F,
    pub eval_0: F,
    pub eval_1: F,
    pub x0: F,
    pub x1: F,
    pub roll_in_value: F,
    pub mmcs_verified: F,
    pub log_height: F,
}

/// Generate trace rows for the CommitPhase table.
///
/// This function takes the FRI folding steps captured from the actual verifier
/// and generates the corresponding AIR trace data.
pub fn generate_commit_phase_trace<F: Field>(steps: Vec<CommitPhaseStep<F>>) -> RowMajorMatrix<F> {
    let num_rows = steps.len();
    let num_cols = core::mem::size_of::<CommitPhaseCols<u8>>();

    let mut trace_data = Vec::with_capacity(num_rows * num_cols);

    for step in steps {
        // Compute intermediate values for the FRI folding formula
        let beta_minus_x0 = step.beta - step.x0;
        let eval_1_minus_eval_0 = step.eval_1 - step.eval_0;
        let x1_minus_x0 = step.x1 - step.x0;

        // Compute inverse: (x1 - x0)^(-1)
        let inverse = x1_minus_x0.inverse();

        // Compute intermediate: (beta - x0) * (eval_1 - eval_0) * (x1 - x0)^(-1)
        let intermediate = beta_minus_x0 * eval_1_minus_eval_0 * inverse;

        // Compute pre-roll-in result: eval_0 + intermediate
        let folded_result_pre = step.eval_0 + intermediate;

        // Final result with roll-in: folded_result_pre + roll_in_value
        let folded_result = folded_result_pre + step.roll_in_value;

        // Create the row data in the same order as CommitPhaseCols struct
        let row_data = [
            step.query_id,
            step.phase_idx,
            step.beta,
            step.domain_index,
            step.parent_index,
            step.sibling_index,
            step.eval_0,
            step.eval_1,
            step.x0,
            step.x1,
            beta_minus_x0,
            eval_1_minus_eval_0,
            x1_minus_x0,
            inverse,
            intermediate,
            folded_result_pre,
            step.roll_in_value,
            folded_result,
            step.mmcs_verified,
            step.log_height,
        ];

        trace_data.extend_from_slice(&row_data);
    }

    RowMajorMatrix::new(trace_data, num_cols)
}

/// Helper function to create a CommitPhaseStep from FRI verifier data.
/// This would be called during actual FRI verification to capture the immediate values.
pub fn create_commit_phase_step<F: Field>(
    query_id: usize,
    phase_idx: usize,
    beta: F,
    domain_index: usize,
    eval_0: F,
    eval_1: F, // sibling value
    x0: F,     // subgroup point
    roll_in_value: F,
    log_height: usize,
    mmcs_verified: bool,
) -> CommitPhaseStep<F> {
    CommitPhaseStep {
        query_id: F::from_usize(query_id),
        phase_idx: F::from_usize(phase_idx),
        beta,
        domain_index: F::from_usize(domain_index),
        parent_index: F::from_usize(domain_index >> 1),
        sibling_index: F::from_usize(domain_index ^ 1),
        eval_0,
        eval_1,
        x0,
        x1: -x0, // x1 = -x0 for arity-2 folding
        roll_in_value,
        mmcs_verified: if mmcs_verified { F::ONE } else { F::ZERO },
        log_height: F::from_usize(log_height),
    }
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
        let step = create_commit_phase_step::<TestField>(
            0,                        // query_id
            0,                        // phase_idx
            TestField::from_u32(42),  // beta
            8,                        // domain_index
            TestField::from_u32(100), // eval_0
            TestField::from_u32(200), // eval_1
            TestField::from_u32(7),   // x0
            TestField::ZERO,          // roll_in_value
            4,                        // log_height
            true,                     // mmcs_verified
        );

        // Verify computed fields
        assert_eq!(step.parent_index, TestField::from_usize(4)); // 8 >> 1
        assert_eq!(step.sibling_index, TestField::from_usize(9)); // 8 ^ 1
        assert_eq!(step.x1, -TestField::from_u32(7)); // -x0
    }

    #[test]
    fn test_trace_generation() {
        let steps = vec![
            create_commit_phase_step::<TestField>(
                0,
                0,
                TestField::from_u32(42),
                8,
                TestField::from_u32(100),
                TestField::from_u32(200),
                TestField::from_u32(7),
                TestField::ZERO,
                4,
                true,
            ),
            create_commit_phase_step::<TestField>(
                0,
                1,
                TestField::from_u32(43),
                4,
                TestField::from_u32(150),
                TestField::from_u32(250),
                TestField::from_u32(3),
                TestField::ZERO,
                3,
                true,
            ),
        ];

        let trace = generate_commit_phase_trace(steps);
        assert_eq!(trace.height(), 2);
        assert_eq!(trace.width(), core::mem::size_of::<CommitPhaseCols<u8>>());
    }
}
