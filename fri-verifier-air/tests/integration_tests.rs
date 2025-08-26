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

    // 1) Generate REAL commit-phase rows from a fresh FRI proof (verifier view)
    let polynomial_log_sizes = [4, 6]; // small for speed
    let commit_phase_steps =
        extract_commit_phase_steps_from_fri_verification(&mut rng, &polynomial_log_sizes);

    // 2) Build the AIR trace from rows
    let mut trace = generate_commit_phase_trace(&commit_phase_steps);

    // 3) Pad to power-of-two height for STARK proving (repeat last row)
    let h = trace.height();
    if h > 0 && !h.is_power_of_two() {
        let target = h.next_power_of_two();
        let pad_rows = target - h;

        let mut padded = trace.values.clone();
        if let Some(last) = trace.row_slice(h - 1) {
            for _ in 0..pad_rows {
                padded.extend_from_slice(&*last);
            }
        }
        trace = RowMajorMatrix::new(padded, trace.width());
    }

    // 4) STARK config (PCS + challenger)
    let config = create_stark_config(&mut rng);

    // 5) Prove & verify with the updated AIR
    let air = CommitPhaseAir::<Val>::new();
    let public_inputs = vec![];

    let proof = prove(&config, &air, trace, &public_inputs);
    verify(&config, &air, &proof, &public_inputs)
        .expect("FRI Verifier AIR STARK proof verification failed");
}
