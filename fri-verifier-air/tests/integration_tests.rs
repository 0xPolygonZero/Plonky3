use p3_matrix::{Matrix, dense::RowMajorMatrix};
use p3_uni_stark::{prove, verify};

use p3_fri_verifier_air::{CommitPhaseAir, generate_commit_phase_trace};

use rand::SeedableRng;
use rand::rngs::SmallRng;

mod proof_generation_utils;
use proof_generation_utils::*;

#[test]
fn test_fri_verifier_air_stark_proving() {
    let mut rng = SmallRng::seed_from_u64(42);

    // Generate real FRI proof data for CommitPhase AIR using utility function
    let polynomial_log_sizes = [4, 6]; // Small sizes for faster testing
    let commit_phase_steps =
        extract_commit_phase_steps_from_fri_verification(&mut rng, &polynomial_log_sizes);

    // Generate trace from real data
    let mut trace = generate_commit_phase_trace(&commit_phase_steps);

    // Pad trace to power of 2 for STARK proving
    let trace_height = trace.height();
    if trace_height > 0 && (trace_height & (trace_height - 1)) != 0 {
        // Not a power of 2, find next power of 2
        let next_power_of_2 = 1 << (trace_height.ilog2() + 1);

        // Pad with dummy rows
        let padding_rows = next_power_of_2 - trace_height;
        let mut padded_data = trace.values.clone();

        // Add dummy rows (just copy the last row)
        if let Some(last_row) = trace.row_slice(trace_height - 1) {
            for _ in 0..padding_rows {
                padded_data.extend_from_slice(&*last_row);
            }
        }

        trace = RowMajorMatrix::new(padded_data, trace.width());
    }

    // Setup STARK configuration using utility function
    let config = create_stark_config(&mut rng);

    // Create AIR instance
    let air = CommitPhaseAir::<Val>::new();

    // Prove
    let proof = prove(&config, &air, trace, &vec![]);

    // Verify
    verify(&config, &air, &proof, &vec![])
        .expect("FRI Verifier AIR STARK proof verification failed");
}
