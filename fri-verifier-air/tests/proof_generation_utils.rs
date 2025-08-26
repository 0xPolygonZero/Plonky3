use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_challenger::{
    CanObserve, CanSampleBits, DuplexChallenger, FieldChallenger,
};
use p3_commit::{ExtensionMmcs, Pcs};
use p3_dft::Radix2Dit;
use p3_field::{Field, PrimeCharacteristicRing, TwoAdicField, extension::BinomialExtensionField, BasedVectorSpace};
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
/// Sets up a PCS, produces a proof, and converts the commit phase info
/// into `CommitPhaseStep` rows compatible with the updated columns/AIR.
pub fn extract_commit_phase_steps_from_fri_verification<R: Rng>(
    rng: &mut R,
    polynomial_log_sizes: &[u8],
) -> Vec<CommitPhaseStep<Val>> {
    // --- PCS setup ---
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

    // Convert polynomial sizes to field elements for transcript
    let val_sizes: Vec<Val> = polynomial_log_sizes
        .iter()
        .map(|&i| Val::from_u8(i))
        .collect();
    let num_evals = polynomial_log_sizes.len();

    // --- Prover world: commit & open ---
    let (commitment, _opened_values, opening_proof) = {
        let mut p_chal = Challenger::new(perm.clone());
        p_chal.observe_slice(&val_sizes);

        // Random evaluation matrices for each degree
        let evaluations: Vec<_> = polynomial_log_sizes
            .iter()
            .map(|deg_bits| {
                let deg = 1 << deg_bits;
                (
                    <MyPcs as Pcs<Challenge, Challenger>>::natural_domain_for_degree(&pcs, deg),
                    // fixed width for convenience
                    RowMajorMatrix::<Val>::rand_nonzero(rng, deg, 16),
                )
            })
            .collect();

        let (commitment, prover_data) =
            <MyPcs as Pcs<Challenge, Challenger>>::commit(&pcs, evaluations);

        // bind commitment then sample opening point
        p_chal.observe(commitment);
        let zeta: Challenge = p_chal.sample_algebra_element();

        // open all polys at zeta
        let open_data = vec![(&prover_data, vec![vec![zeta]; num_evals])];
        let (opened_values, opening_proof) = pcs.open(open_data, &mut p_chal);

        (commitment, opened_values, opening_proof)
    };

    // --- Turn proof into CommitPhase steps (verifier view) ---
    extract_commit_phase_steps_from_opening_proof(
        &opening_proof,
        &pcs,
        &perm,
        &val_sizes,
        commitment,
    )
}

/// Extract CommitPhase steps from an existing FRI opening proof (verifier-side).
///
/// This mirrors verifier transcript ops:
///  - observe sizes, commitment
///  - sample zeta, alpha
///  - observe commit_phase_commits and sample betas
///  - observe final poly
///  - sample random index
/// Then constructs per-phase folding rows with zero roll-ins (for this util).
fn extract_commit_phase_steps_from_opening_proof(
    opening_proof: &<MyPcs as Pcs<Challenge, Challenger>>::Proof,
    pcs: &MyPcs,
    perm: &Perm,
    val_sizes: &[Val],
    commitment: <MyPcs as Pcs<Challenge, Challenger>>::Commitment,
) -> Vec<CommitPhaseStep<Val>> {
    let mut v_chal = Challenger::new(perm.clone());

    // Bind sizes & commitment; sample opening point
    v_chal.observe_slice(val_sizes);
    v_chal.observe(commitment);
    let _zeta: Challenge = v_chal.sample_algebra_element();

    // Alpha (batch combiner) — not used directly in commit-phase rows, but sampled by verifier
    let _alpha: Challenge = v_chal.sample_algebra_element();

    // Betas: observe each commit, then sample βᵢ
    let betas: Vec<Challenge> = opening_proof
        .commit_phase_commits
        .iter()
        .map(|comm| {
            v_chal.observe(comm.clone());
            v_chal.sample_algebra_element()
        })
        .collect();

    // Bind final polynomial coefficients
    for coeff in &opening_proof.final_poly {
        v_chal.observe_algebra_element(*coeff);
    }

    // Random query index with log_max_height bits
    let log_max_height = opening_proof.commit_phase_commits.len()
        + pcs.fri_params().log_blowup
        + pcs.fri_params().log_final_poly_len;
    let mut domain_index = v_chal.sample_bits(log_max_height);

    // We use the first query proof (single-query typical)
    let query_proof = &opening_proof.query_proofs[0];

    // Start folded_eval at 0 (no top-level reduced opening provided in this util)
    let mut current_folded = Challenge::ZERO;

    let mut steps = Vec::with_capacity(betas.len());
    for (phase_idx, (commit, opening)) in opening_proof
        .commit_phase_commits
        .iter()
        .zip(query_proof.commit_phase_openings.iter())
        .enumerate()
    {
        let beta = betas[phase_idx];
        let sibling = opening.sibling_value;

        // Height and subgroup points for this phase (arity=2 folding)
        let log_folded_height = log_max_height - phase_idx - 1;
        let rev_bits = p3_util::reverse_bits_len(domain_index, log_folded_height);
        let generator = Val::two_adic_generator(log_folded_height + 1);
        let x0 = generator.exp_u64(rev_bits as u64);
        let x1 = -x0;

        // eval_0 is current folded value; eval_1 is sibling from proof
        let eval_0 = current_folded;
        let eval_1 = sibling;

        // β - x0
        let x0_ext =
            Challenge::from_basis_coefficients_slice(&[x0, Val::ZERO, Val::ZERO, Val::ZERO])
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
        let x_diff_inv_ext = Challenge::from_basis_coefficients_slice(&[
            x_diff_inv,
            Val::ZERO,
            Val::ZERO,
            Val::ZERO,
        ])
        .expect("valid inv ext");
        let interpolation_term = beta_eval_product * x_diff_inv_ext;

        // Folded eval BEFORE roll-in
        let folded_eval_pre_rollin = eval_0 + interpolation_term;

        // This util sets roll-ins to zero (is_roll_in=0); you can wire real reduced openings here.
        let beta_squared = beta * beta;
        let is_roll_in = Val::ZERO;
        let reduced_opening = Challenge::ZERO;
        let roll_in_contribution = Challenge::ZERO;

        // Final folded value for this phase (no roll-in)
        let folded_eval = folded_eval_pre_rollin + roll_in_contribution;

        // Commit digest for the phase
        let fri_commit: [Val; 8] = *commit.as_ref();

        // Push row
        steps.push(CommitPhaseStep {
            query_index: Val::from_usize(0),
            phase_index: Val::from_usize(phase_idx),

            beta,
            eval_0,
            eval_1,

            is_roll_in,
            reduced_opening,
            beta_squared,
            roll_in_contribution,

            domain_index: Val::from_usize(domain_index),
            parent_index: Val::from_usize(domain_index >> 1),
            sibling_index: Val::from_usize(domain_index ^ 1),
            lsb: Val::from_usize(domain_index & 1),

            log_height: Val::from_usize(log_folded_height),

            x0,
            x1,

            beta_minus_x0,
            eval_diff,
            x_diff,
            x_diff_inv,
            beta_eval_product,
            interpolation_term,
            folded_eval,

            fri_commit,
        });

        // Move to parent index for next fold
        domain_index >>= 1;
        // Update carry value
        current_folded = folded_eval;
    }

    steps
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
