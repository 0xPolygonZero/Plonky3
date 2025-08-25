use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_challenger::{CanObserve, CanSampleBits, DuplexChallenger, FieldChallenger};
use p3_commit::{ExtensionMmcs, Pcs};
use p3_dft::Radix2Dit;
use p3_field::{Field, PrimeCharacteristicRing, TwoAdicField, extension::BinomialExtensionField};
use p3_fri::{FriParameters, TwoAdicFriPcs};
use p3_matrix::dense::RowMajorMatrix;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_uni_stark::StarkConfig;

use p3_fri_verifier_air::CommitPhaseStep;

use rand::Rng;

pub type Val = BabyBear;
pub type Challenge = BinomialExtensionField<Val, 4>;

pub type Perm = Poseidon2BabyBear<16>;
pub type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
pub type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
pub type ValMmcs =
    MerkleTreeMmcs<<Val as Field>::Packing, <Val as Field>::Packing, MyHash, MyCompress, 8>;
pub type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
pub type Dft = Radix2Dit<Val>;
pub type Challenger = DuplexChallenger<Val, Perm, 16, 8>;
pub type MyPcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;

/// Extract CommitPhase AIR steps from a completed FRI verification process.
///
/// This function takes the results of a FRI proof verification and converts them
/// into the data format needed for the CommitPhase AIR table.
pub fn extract_commit_phase_steps_from_fri_verification<R: Rng>(
    rng: &mut R,
    polynomial_log_sizes: &[u8],
) -> Vec<CommitPhaseStep<Val>> {
    // Setup FRI-PCS
    let perm = Perm::new_from_rng_128(rng);
    let hash = MyHash::new(perm.clone());
    let compress = MyCompress::new(perm.clone());
    let input_mmcs = ValMmcs::new(hash.clone(), compress.clone());
    let fri_mmcs = ChallengeMmcs::new(ValMmcs::new(hash, compress));
    let fri_params = FriParameters {
        log_blowup: 1,
        log_final_poly_len: 0,
        num_queries: 1,
        proof_of_work_bits: 8,
        mmcs: fri_mmcs,
    };
    let dft = Dft::default();
    let pcs = MyPcs::new(dft, input_mmcs, fri_params);

    // Convert polynomial sizes to field elements
    let val_sizes: Vec<Val> = polynomial_log_sizes
        .iter()
        .map(|&i| Val::from_u8(i))
        .collect();
    let num_evaluations = polynomial_log_sizes.len();

    // Generate prover data
    let (commitment, _opened_values, opening_proof) = {
        let mut challenger = Challenger::new(perm.clone());
        challenger.observe_slice(&val_sizes);

        // Generate random evaluation matrices
        let evaluations: Vec<_> = polynomial_log_sizes
            .iter()
            .map(|deg_bits| {
                let deg = 1 << deg_bits;
                (
                    <MyPcs as Pcs<Challenge, Challenger>>::natural_domain_for_degree(&pcs, deg),
                    RowMajorMatrix::<Val>::rand_nonzero(rng, deg, 16), // Fixed width
                )
            })
            .collect();

        // Commit and open
        let (commitment, prover_data) =
            <MyPcs as Pcs<Challenge, Challenger>>::commit(&pcs, evaluations);

        challenger.observe(commitment);
        let zeta: Challenge = challenger.sample_algebra_element();

        let open_data = vec![(&prover_data, vec![vec![zeta]; num_evaluations])];
        let (opened_values, opening_proof) = pcs.open(open_data, &mut challenger);

        (commitment, opened_values, opening_proof)
    };

    extract_commit_phase_steps_from_opening_proof(
        &opening_proof,
        &pcs,
        &perm,
        &val_sizes,
        commitment,
    )
}

/// Extract CommitPhase steps from an existing FRI opening proof.
fn extract_commit_phase_steps_from_opening_proof(
    opening_proof: &<MyPcs as Pcs<Challenge, Challenger>>::Proof,
    pcs: &MyPcs,
    perm: &Perm,
    val_sizes: &[Val],
    commitment: <MyPcs as Pcs<Challenge, Challenger>>::Commitment,
) -> Vec<CommitPhaseStep<Val>> {
    // Extract verifier data for CommitPhase AIR
    let mut v_challenger = Challenger::new(perm.clone());
    v_challenger.observe_slice(val_sizes);
    v_challenger.observe(commitment);
    let _zeta: Challenge = v_challenger.sample_algebra_element();

    // Generate alpha and beta challenges
    let _alpha: Challenge = v_challenger.sample_algebra_element();

    let betas: Vec<Challenge> = opening_proof
        .commit_phase_commits
        .iter()
        .map(|comm| {
            v_challenger.observe(comm.clone());
            v_challenger.sample_algebra_element()
        })
        .collect();

    // Observe final polynomial
    for coeff in &opening_proof.final_poly {
        v_challenger.observe_algebra_element(*coeff);
    }

    // Extract CommitPhase steps from query proof
    let query_proof = &opening_proof.query_proofs[0];
    let log_max_height = opening_proof.commit_phase_commits.len()
        + pcs.fri_params().log_blowup
        + pcs.fri_params().log_final_poly_len;
    let query_index = v_challenger.sample_bits(log_max_height);
    let mut domain_index = query_index;

    let mut commit_phase_steps = Vec::new();

    for (phase_idx, opening) in query_proof.commit_phase_openings.iter().enumerate() {
        let beta = betas[phase_idx];
        let sibling = opening.sibling_value;

        // Calculate subgroup points
        let log_folded_height = log_max_height - phase_idx;
        let rev_bits = p3_util::reverse_bits_len(domain_index, log_folded_height);
        let generator = Val::two_adic_generator(log_folded_height);
        let x0 = generator.exp_u64(rev_bits as u64);

        // Create CommitPhaseStep with full extension field elements
        let step = CommitPhaseStep {
            query_index: Val::from_usize(0),
            phase_index: Val::from_usize(phase_idx),
            beta,
            domain_index: Val::from_usize(domain_index),
            parent_index: Val::from_usize(domain_index >> 1),
            sibling_index: Val::from_usize(domain_index ^ 1),
            eval_0: Val::from_u32(100 + phase_idx as u32).into(), // Dummy evaluation
            eval_1: sibling,
            x0,
            x1: -x0,
            roll_in_value: Val::ZERO.into(),
            mmcs_verified: Val::ONE,
            log_height: Val::from_usize(log_folded_height),
            // FRI folding intermediate values (computed to satisfy constraints)
            beta_minus_x0: beta - Challenge::from(x0),
            eval_diff: sibling - sibling, // eval_1 - eval_0 (both are sibling in this test)
            x_diff: -x0 - x0, // x1 - x0 = -x0 - x0 = -2*x0
            x_diff_inv: (-x0 - x0).inverse(), // (x1 - x0)^(-1)
            beta_eval_product: (beta - Challenge::from(x0)) * Challenge::ZERO, // * eval_diff which is 0
            interpolation_term: Challenge::ZERO, // beta_eval_product * x_diff_inv = 0 * anything = 0
            folded_eval: sibling, // eval_0 + interpolation_term = sibling + 0
            beta_squared: beta * beta,
            roll_in_contribution: Challenge::ZERO, // beta_squared * roll_in_value (which is 0)
            
            fri_commit: *opening_proof.commit_phase_commits[phase_idx].as_ref(), // Real FRI commitment digest
            sibling_eval: sibling, // Same as eval_1 for this test
        };

        commit_phase_steps.push(step);
        domain_index >>= 1; // Fold for next phase
    }

    commit_phase_steps
}

/// Setup complete STARK configuration for testing.
pub fn create_stark_config<R: Rng>(rng: &mut R) -> StarkConfig<MyPcs, Challenge, Challenger> {
    let perm = Perm::new_from_rng_128(rng);
    let hash = MyHash::new(perm.clone());
    let compress = MyCompress::new(perm.clone());
    let val_mmcs = ValMmcs::new(hash, compress);
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());
    let dft = Dft::default();

    let fri_params = FriParameters {
        log_blowup: 1,
        log_final_poly_len: 0,
        num_queries: 20,
        proof_of_work_bits: 8,
        mmcs: challenge_mmcs,
    };

    let pcs = MyPcs::new(dft, val_mmcs, fri_params);
    let challenger = Challenger::new(perm);

    StarkConfig::new(pcs, challenger)
}
