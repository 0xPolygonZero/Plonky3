use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_challenger::{
    CanObserve, CanSampleBits, DuplexChallenger, FieldChallenger, GrindingChallenger,
};
use p3_commit::{ExtensionMmcs, Mmcs, Pcs};
use p3_dft::Radix2Dit;
use p3_field::coset::TwoAdicMultiplicativeCoset;
use p3_field::extension::BinomialExtensionField;
use p3_field::{Field, PrimeCharacteristicRing, TwoAdicField};
use p3_fri::{FriParameters, TwoAdicFriPcs};
use p3_matrix::dense::RowMajorMatrix;
use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

type Val = BabyBear;
type Challenge = BinomialExtensionField<Val, 4>;

type Perm = Poseidon2BabyBear<16>;
type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
type ValMmcs =
    MerkleTreeMmcs<<Val as Field>::Packing, <Val as Field>::Packing, MyHash, MyCompress, 8>;
type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
type Challenger = DuplexChallenger<Val, Perm, 16, 8>;
type MyPcs = TwoAdicFriPcs<BabyBear, Radix2Dit<BabyBear>, ValMmcs, ChallengeMmcs>;

/// Returns a permutation and a FRI-pcs instance.
fn get_ldt_for_testing<R: Rng>(rng: &mut R, log_final_poly_len: usize) -> (Perm, MyPcs) {
    let perm = Perm::new_from_rng_128(rng);
    let hash = MyHash::new(perm.clone());
    let compress = MyCompress::new(perm.clone());
    let input_mmcs = ValMmcs::new(hash.clone(), compress.clone());
    let fri_mmcs = ChallengeMmcs::new(ValMmcs::new(hash, compress));
    let fri_params = FriParameters {
        log_blowup: 1,
        log_final_poly_len,
        num_queries: 1,
        proof_of_work_bits: 8,
        mmcs: fri_mmcs,
    };
    let dft = Radix2Dit::default();
    let pcs = MyPcs::new(dft, input_mmcs, fri_params);
    (perm, pcs)
}

/// Check that the loop of `pcs.commit`, `pcs.open`, and `pcs.verify` work correctly.
///
/// We create a random polynomial of size `1 << log_size` for each size in `polynomial_log_sizes`.
/// We then commit to these polynomials using a `log_blowup` of `1`.
///
/// We open each polynomial at the same point `zeta` and run FRI to verify the openings, stopping
/// FRI at `log_final_poly_len`.
fn do_test_fri_ldt<R: Rng>(rng: &mut R, log_final_poly_len: usize, polynomial_log_sizes: &[u8]) {
    let (perm, pcs) = get_ldt_for_testing(rng, log_final_poly_len);

    // Convert the polynomial_log_sizes into field elements so they can be observed.
    let val_sizes: Vec<Val> = polynomial_log_sizes
        .iter()
        .map(|&i| Val::from_u8(i))
        .collect();

    // --- Prover World ---
    let (commitment, opened_values, opening_proof, mut p_challenger) = {
        // Initialize the challenger and observe the `polynomial_log_sizes`.
        let mut challenger = Challenger::new(perm.clone());
        challenger.observe_slice(&val_sizes);

        // Generate random evaluation matrices for each polynomial degree.
        let evaluations: Vec<(TwoAdicMultiplicativeCoset<Val>, RowMajorMatrix<Val>)> =
            polynomial_log_sizes
                .iter()
                .map(|deg_bits| {
                    let deg = 1 << deg_bits;
                    (
                        // Get the TwoAdicSubgroup of this degree.
                        <MyPcs as Pcs<Challenge, Challenger>>::natural_domain_for_degree(&pcs, deg),
                        // Generate a random matrix of evaluations.
                        RowMajorMatrix::<Val>::rand_nonzero(rng, deg, 16),
                    )
                })
                .collect();

        let num_evaluations = evaluations.len();

        // Commit to all the evaluation matrices.
        let (commitment, prover_data) =
            <TwoAdicFriPcs<BabyBear, Radix2Dit<BabyBear>, ValMmcs, ChallengeMmcs> as Pcs<
                Challenge,
                Challenger,
            >>::commit(&pcs, evaluations);

        // Observe the commitment.
        challenger.observe(commitment);

        // Sample the challenge point zeta which all polynomials
        // will be opened at.
        let zeta: Challenge = challenger.sample_algebra_element();

        // Prepare the data into the form expected by `pcs.open`.
        let open_data = vec![(&prover_data, vec![vec![zeta]; num_evaluations])]; // open every chunk at zeta

        // Open all polynomials at zeta and produce the opening proof.
        let (opened_values, opening_proof) = pcs.open(open_data, &mut challenger);

        // Return the commitment, opened values, opening proof and challenger.
        // The first three of these are always passed to the verifier. The
        // last is to double check that the prover and verifiers challengers
        // agree at appropriate points.
        (commitment, opened_values, opening_proof, challenger)
    };

    // --- Verifier World ---
    let mut v_challenger = {
        // Initialize the verifier's challenger with the same permutation.
        // Observe the `polynomial_log_sizes` and `commitment` in the same order
        // as the prover.
        let mut challenger = Challenger::new(perm.clone());
        challenger.observe_slice(&val_sizes);
        challenger.observe(commitment);

        // Sample the opening point.
        let zeta = challenger.sample_algebra_element();

        // Construct the expected initial polynomial domains.
        // Right now it doesn't matter what these are so long as the size
        // is correct.
        let domains = polynomial_log_sizes.iter().map(|&size| {
            <MyPcs as Pcs<Challenge, Challenger>>::natural_domain_for_degree(&pcs, 1 << size)
        });

        // Prepare the data into the form expected by `pcs.verify`.
        // Note that commitment and opened_values are always sent by
        // the prover.
        let commitments_with_opening_points = vec![(
            commitment,
            domains
                .into_iter()
                .zip(opened_values.into_iter().flatten().flatten())
                .map(|(domain, value)| (domain, vec![(zeta, value)]))
                .collect(),
        )];

        // Verify the opening proof.
        let verification = pcs.verify(
            commitments_with_opening_points,
            &opening_proof,
            &mut challenger,
        );
        assert!(verification.is_ok());
        challenger
    };

    // Check that the prover and verifier challengers agree.
    assert_eq!(
        p_challenger.sample_bits(8),
        v_challenger.sample_bits(8),
        "prover and verifier transcript have same state after FRI"
    );
}

/// Test that the FRI commit, open and verify process work correctly
/// for a range of `final_poly_degree` values.
#[test]
fn test_fri_ldt() {
    // Chosen to ensure there are both multiple polynomials
    // of the same size and that the array is not ordered.
    let polynomial_log_sizes = [5, 8, 10, 7, 5, 5, 7];
    for i in 0..5 {
        let mut rng = SmallRng::seed_from_u64(i as u64);
        do_test_fri_ldt(&mut rng, i, &polynomial_log_sizes);
    }
}

/// This test is expected to panic because there is a polynomial degree which
/// the prover commits too which is less than `final_poly_degree`.
#[test]
#[should_panic]
fn test_fri_ldt_should_panic() {
    // Chosen to ensure there are both multiple polynomials
    // of the same size and that the array is not ordered.
    let polynomial_log_sizes = [5, 8, 10, 7, 5, 5, 7];
    let mut rng = SmallRng::seed_from_u64(5);
    do_test_fri_ldt(&mut rng, 5, &polynomial_log_sizes);
}

/// Test that outputs all immediate values from FRI verifier for AIR table design
#[test]
fn test_fri_verifier_immediate_values() {
    let polynomial_log_sizes = [5, 8, 8, 10];
    let mut rng = SmallRng::seed_from_u64(0);

    // Create a custom test function that captures all immediate values
    do_test_fri_verifier_immediate_values(&mut rng, 0, &polynomial_log_sizes);
}

/// Test function that captures all immediate values from FRI verifier
fn do_test_fri_verifier_immediate_values<R: Rng>(
    rng: &mut R,
    log_final_poly_len: usize,
    polynomial_log_sizes: &[u8],
) {
    let (perm, pcs) = get_ldt_for_testing(rng, log_final_poly_len);

    // Convert the polynomial_log_sizes into field elements so they can be observed.
    let val_sizes: Vec<Val> = polynomial_log_sizes
        .iter()
        .map(|&i| Val::from_u8(i))
        .collect();

    // --- Prover World ---
    let (commitment, opened_values, opening_proof, mut p_challenger) = {
        // Initialize the challenger and observe the `polynomial_log_sizes`.
        let mut challenger = Challenger::new(perm.clone());
        challenger.observe_slice(&val_sizes);

        // Generate random evaluation matrices for each polynomial degree.
        let evaluations: Vec<(TwoAdicMultiplicativeCoset<Val>, RowMajorMatrix<Val>)> =
            polynomial_log_sizes
                .iter()
                .map(|deg_bits| {
                    let deg = 1 << deg_bits;
                    (
                        // Get the TwoAdicSubgroup of this degree.
                        <MyPcs as Pcs<Challenge, Challenger>>::natural_domain_for_degree(&pcs, deg),
                        // Generate a random matrix of evaluations.
                        RowMajorMatrix::<Val>::rand_nonzero(rng, deg, *deg_bits as usize - 4),
                    )
                })
                .collect();

        let num_evaluations = evaluations.len();

        // Commit to all the evaluation matrices.
        let (commitment, prover_data) =
            <TwoAdicFriPcs<BabyBear, Radix2Dit<BabyBear>, ValMmcs, ChallengeMmcs> as Pcs<
                Challenge,
                Challenger,
            >>::commit(&pcs, evaluations);

        // Observe the commitment.
        challenger.observe(commitment);

        // Sample the challenge point zeta which all polynomials
        // will be opened at.
        let zeta: Challenge = challenger.sample_algebra_element();

        // Prepare the data into the form expected by `pcs.open`.
        let open_data = vec![(&prover_data, vec![vec![zeta]; num_evaluations])]; // open every chunk at zeta

        // Open all polynomials at zeta and produce the opening proof.
        let (opened_values, opening_proof) = pcs.open(open_data, &mut challenger);

        // Return the commitment, opened values, opening proof and challenger.
        // The first three of these are always passed to the verifier. The
        // last is to double check that the prover and verifiers challengers
        // agree at appropriate points.
        (commitment, opened_values, opening_proof, challenger)
    };

    // --- Verifier World ---
    let mut v_challenger = {
        // Initialize the verifier's challenger with the same permutation.
        // Observe the `polynomial_log_sizes` and `commitment` in the same order
        // as the prover.
        let mut challenger = Challenger::new(perm.clone());
        challenger.observe_slice(&val_sizes);
        challenger.observe(commitment);

        // Sample the opening point.
        let zeta = challenger.sample_algebra_element();

        // Construct the expected initial polynomial domains.
        // Right now it doesn't matter what these are so long as the size
        // is correct.
        let domains = polynomial_log_sizes.iter().map(|&size| {
            <MyPcs as Pcs<Challenge, Challenger>>::natural_domain_for_degree(&pcs, 1 << size)
        });
        println!("Domains: {:?}", domains);

        // Prepare the data into the form expected by `pcs.verify`.
        // Note that commitment and opened_values are always sent by
        // the prover.
        let commitments_with_opening_points = vec![(
            commitment,
            domains
                .into_iter()
                .zip(opened_values.into_iter().flatten().flatten())
                .map(|(domain, value)| (domain, vec![(zeta, value)]))
                .collect(),
        )];

        assert!(
            capture_fri_verifier_immediate_values(
                commitments_with_opening_points.clone(),
                &opening_proof,
                &mut challenger,
                &pcs, // Pass the PCS instance for real verification calls
            )
            .is_ok()
        );

        challenger
    };

    // Check that the prover and verifier challengers agree.
    assert_eq!(
        p_challenger.sample_bits(8),
        v_challenger.sample_bits(8),
        "prover and verifier transcript have same state after FRI"
    );
}

/// Function to capture all immediate values from FRI verifier with REAL verification calls
fn capture_fri_verifier_immediate_values(
    commitments_with_opening_points: Vec<(
        <MyPcs as Pcs<Challenge, Challenger>>::Commitment,
        Vec<(
            TwoAdicMultiplicativeCoset<Val>,
            Vec<(Challenge, Vec<Challenge>)>,
        )>,
    )>,
    opening_proof: &<MyPcs as Pcs<Challenge, Challenger>>::Proof,
    challenger: &mut Challenger,
    pcs: &MyPcs, // Added to access MMCS instances and parameters
) -> Result<(), Box<dyn std::error::Error>> {
    // CRITICAL: First step in PCS verify - observe all evaluation values
    // This happens BEFORE calling verify_fri in the real PCS implementation
    println!("=== PCS VERIFIER STEP: OBSERVING EVALUATIONS ===");
    for (_, round) in &commitments_with_opening_points {
        for (_, mat) in round {
            for (_, point) in mat {
                point.iter().for_each(|&opening| {
                    challenger.observe_algebra_element(opening);
                    println!("  Observed evaluation: {:?}", opening);
                });
            }
        }
    }

    // Extract FRI proof from opening proof
    let _fri_proof = match opening_proof {
        p3_fri::FriProof {
            commit_phase_commits,
            query_proofs,
            final_poly,
            pow_witness,
        } => {
            println!("\n=== FRI VERIFIER IMMEDIATE VALUES ===");
            println!("=== FRI PARAMETERS ===");
            println!(
                "Number of commit phase commitments: {}",
                commit_phase_commits.len()
            );
            println!("Number of query proofs: {}", query_proofs.len());
            println!("Final polynomial degree: {}", final_poly.len() - 1);
            println!("Proof of work witness: {:?}", pow_witness);

            println!("\n=== COMMIT PHASE COMMITMENTS ===");
            println!(
                "Number of commit phase commitments: {}",
                commit_phase_commits.len()
            );
            for (i, commit) in commit_phase_commits.iter().enumerate() {
                println!("Commit phase {}: {:?}", i, commit);
            }

            println!("\n=== FINAL POLYNOMIAL ===");
            println!("Final polynomial coefficients: {:?}", final_poly);
            println!("Final polynomial degree: {}", final_poly.len() - 1);

            println!("\n=== PROOF OF WORK ===");
            println!("PoW witness: {:?}", pow_witness);

            println!("\n=== QUERY PROOFS ===");
            println!("Number of query proofs: {}", query_proofs.len());

            for (query_idx, query_proof) in query_proofs.iter().enumerate() {
                println!("\n--- Query Proof {} ---", query_idx);

                println!("Input proof length: {}", query_proof.input_proof.len());
                println!(
                    "Number of commit phase openings: {}",
                    query_proof.commit_phase_openings.len()
                );

                for (phase_idx, opening) in query_proof.commit_phase_openings.iter().enumerate() {
                    println!("  Phase {} opening:", phase_idx);
                    println!("    Sibling value: {:?}", opening.sibling_value);
                    println!(
                        "    Opening proof type: {:?}",
                        std::any::type_name_of_val(&opening.opening_proof)
                    );
                }
            }

            // Simulate the verification process to capture more immediate values
            simulate_fri_verification_immediate_values(
                &commit_phase_commits,
                &query_proofs,
                &final_poly,
                &pow_witness,
                challenger,
                commitments_with_opening_points.as_slice(),
                pcs, // Pass the PCS instance
            )?;
        }
    };

    Ok(())
}

/// Simulate FRI verification to capture immediate values - must match exact challenger operations
fn simulate_fri_verification_immediate_values(
    commit_phase_commits: &[<ChallengeMmcs as p3_commit::Mmcs<Challenge>>::Commitment],
    query_proofs: &[p3_fri::QueryProof<
        Challenge,
        ChallengeMmcs,
        Vec<p3_commit::BatchOpening<Val, ValMmcs>>,
    >],
    final_poly: &[Challenge],
    pow_witness: &<Challenger as p3_challenger::GrindingChallenger>::Witness,
    challenger: &mut Challenger,
    _commitments_with_opening_points: &[(
        <MyPcs as Pcs<Challenge, Challenger>>::Commitment,
        Vec<(
            TwoAdicMultiplicativeCoset<Val>,
            Vec<(Challenge, Vec<Challenge>)>,
        )>,
    )],
    pcs: &MyPcs, // Added PCS parameter
) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n=== VERIFICATION SIMULATION IMMEDIATE VALUES ===");

    // Step 1: Generate the Batch combination challenge (line 73 in verifier.rs)
    let alpha: Challenge = challenger.sample_algebra_element();
    println!("Alpha (batch combination challenge): {:?}", alpha);

    // Use REAL FRI parameters from the PCS
    let fri_params = pcs.fri_params();
    let log_blowup = fri_params.log_blowup;
    let log_final_poly_len = fri_params.log_final_poly_len;
    let num_queries = fri_params.num_queries;
    let proof_of_work_bits = fri_params.proof_of_work_bits;
    let log_global_max_height = commit_phase_commits.len() + log_blowup + log_final_poly_len;
    let log_max_height = commit_phase_commits.len() + log_blowup + log_final_poly_len;
    let log_final_height = log_blowup + log_final_poly_len;

    println!("log_global_max_height: {}", log_global_max_height);
    println!("log_max_height: {}", log_max_height);
    println!("log_final_height: {}", log_final_height);

    // Step 2: Generate all of the random challenges for the FRI rounds (lines 82-91)
    // CRITICAL: Must observe commitment BEFORE sampling beta (line 88-89)
    let betas: Vec<Challenge> = commit_phase_commits
        .iter()
        .map(|comm| {
            // To match with the prover (and for security purposes),
            // we observe the commitment before sampling the challenge.
            challenger.observe(comm.clone());
            challenger.sample_algebra_element()
        })
        .collect();

    println!("\n=== BETA CHALLENGES ===");
    for (i, beta) in betas.iter().enumerate() {
        println!("Beta {}: {:?}", i, beta);
    }

    // Step 3: Observe all coefficients of the final polynomial (lines 99-102)
    // VALIDATION: Ensure final polynomial has expected degree (lines 94-96)
    println!("\n=== FINAL POLYNOMIAL VALIDATION ===");
    let expected_final_poly_len = fri_params.final_poly_len(); // Use real FRI method
    if final_poly.len() != expected_final_poly_len {
        return Err(format!(
            "InvalidProofShape: Final polynomial length {} != expected {}",
            final_poly.len(),
            expected_final_poly_len
        )
        .into());
    }
    println!(
        "✓ Final polynomial has expected degree: {} coefficients (constant polynomial)",
        final_poly.len()
    );
    assert_eq!(
        final_poly.len(),
        1,
        "With log_final_poly_len=0, final polynomial should be constant"
    );

    println!("\n=== FINAL POLYNOMIAL OBSERVATION ===");
    println!("Observing final polynomial coefficient (constant):");
    let final_constant = final_poly[0];
    challenger.observe_algebra_element(final_constant);
    println!("  Final polynomial constant: {:?}", final_constant);

    // VALIDATION: Ensure we have expected number of query proofs (lines 105-107)
    if query_proofs.len() != num_queries {
        return Err(format!(
            "InvalidProofShape: Query proofs {} != expected {}",
            query_proofs.len(),
            num_queries
        )
        .into());
    }
    println!("✓ Correct number of query proofs: {}", query_proofs.len());

    // Step 4: REAL Check PoW (Proof of Work) (lines 110-112)
    println!("\n=== PROOF OF WORK CHECK ===");
    println!("PoW bits: {}", proof_of_work_bits);
    println!("PoW witness: {:?}", pow_witness);

    // REAL proof of work verification
    if !challenger.check_witness(proof_of_work_bits, pow_witness.clone()) {
        return Err("InvalidPowWitness: Proof of work verification failed".into());
    }
    println!("✓ Proof of work verification PASSED");

    // Step 5: Process each query proof (lines 121-188)
    for (query_idx, query_proof) in query_proofs.iter().enumerate() {
        println!("\n--- Query {} Verification ---", query_idx);

        // Step 5a: Generate the random index (line 127)
        // Note: extra_query_index_bits() is 0 for standard folding
        let index = challenger.sample_bits(log_max_height);
        println!("Query index: {} ({} bits)", index, log_max_height);

        // Step 5b: Implement REAL input opening verification (open_input function)
        println!("=== INPUT OPENING VERIFICATION (REAL IMPLEMENTATION) ===");

        // Implement the actual open_input logic from verifier.rs:349-442
        // This processes the real input_proof and commitments_with_opening_points
        let mut reduced_openings_map =
            std::collections::BTreeMap::<usize, (Challenge, Challenge)>::new();

        // Process each batch commitment and opening proof
        for (batch_opening, (batch_commit, mats)) in query_proof
            .input_proof
            .iter()
            .zip(_commitments_with_opening_points.iter())
        {
            println!("Processing batch commitment: {:?}", batch_commit);

            // Find the height of each matrix in the batch
            let batch_heights: Vec<usize> = mats
                .iter()
                .map(|(domain, _)| domain.size() << log_blowup)
                .collect();

            println!("Batch heights: {:?}", batch_heights);

            // CRITICAL: REAL MMCS verification of input batch opening (lines 381-388 in verifier.rs)
            let batch_dims: Vec<p3_matrix::Dimensions> = batch_heights
                .iter()
                .map(|&height| p3_matrix::Dimensions { width: 0, height }) // width is 0 for MMCS
                .collect();

            let reduced_index = batch_heights
                .iter()
                .max()
                .map(|&h| index >> (log_global_max_height - h.ilog2() as usize))
                .unwrap_or(0);

            println!(
                "  MMCS verification - batch_dims: {:?}, reduced_index: {}",
                batch_dims, reduced_index
            );

            // REAL MMCS verification call using the actual input MMCS from PCS
            println!("  === MMCS INPUT VERIFICATION IMMEDIATE VALUES ===");
            println!("  Function: input_mmcs.verify_batch()");
            println!("  Inputs:");
            println!("    batch_commit: {:?}", batch_commit);
            println!("    batch_dims: {:?}", batch_dims);
            println!("    reduced_index: {}", reduced_index);
            println!(
                "    batch_opening data: {} matrices",
                batch_opening.opened_values.len()
            );
            for (mat_idx, mat_data) in batch_opening.opened_values.iter().enumerate() {
                println!("      Matrix {}: {} opened values", mat_idx, mat_data.len());
                for (val_idx, val) in mat_data.iter().enumerate() {
                    println!("        Value[{}]: {:?}", val_idx, val);
                }
            }
            println!("  Equation: verify_batch(C, dims, idx, π) → bool");
            println!("    where C = Merkle commitment, dims = matrix dimensions");
            println!("          idx = leaf index, π = opening proof");

            // ACTUAL MMCS verification call with the REAL input MMCS
            let verification_result = pcs.input_mmcs().verify_batch(
                batch_commit,
                &batch_dims,
                reduced_index,
                batch_opening.into(),
            );

            match verification_result {
                Ok(()) => {
                    println!("  Output: SUCCESS (verification passed)");
                    println!("  ✓ MMCS batch verification PASSED - input proof verified");
                }
                Err(e) => {
                    println!("  Output: FAILED - {:?}", e);
                    return Err(format!("InputError: MMCS verification failed: {:?}", e).into());
                }
            }

            // For each matrix in the commitment
            for (mat_opening, (mat_domain, mat_points_and_values)) in
                batch_opening.opened_values.iter().zip(mats.iter())
            {
                let log_height = (mat_domain.size() << log_blowup).ilog2() as usize;

                let bits_reduced = log_global_max_height - log_height;
                let rev_reduced_index =
                    p3_util::reverse_bits_len(index >> bits_reduced, log_height);

                // Compute evaluation point x = g * h^rev_reduced_index
                let x = Val::GENERATOR
                    * Val::two_adic_generator(log_height).exp_u64(rev_reduced_index as u64);

                println!(
                    "Matrix domain size: {}, log_height: {}, x: {:?}",
                    mat_domain.size(),
                    log_height,
                    x
                );

                let (alpha_pow, ro) = reduced_openings_map
                    .entry(log_height)
                    .or_insert((Challenge::ONE, Challenge::ZERO));

                // For each polynomial in the matrix, compute (f(z) - f(x))/(z - x)
                for (z, ps_at_z) in mat_points_and_values {
                    println!("  === REDUCED OPENING COMPUTATION ===");
                    println!("    Evaluation point z: {:?}", z);
                    println!("    Evaluation point x: {:?}", x);
                    println!("    Equation: quotient = (z - x)^(-1)");
                    let quotient = (*z - x).inverse();
                    println!("    quotient = ({:?} - {:?})^(-1) = {:?}", z, x, quotient);

                    for (poly_idx, (&p_at_x, &p_at_z)) in
                        mat_opening.iter().zip(ps_at_z.iter()).enumerate()
                    {
                        println!("    Polynomial {} at height {}:", poly_idx, log_height);
                        println!("      p_at_x (at domain point): {:?}", p_at_x);
                        println!("      p_at_z (at challenge point): {:?}", p_at_z);

                        let contribution = *alpha_pow * (p_at_z - p_at_x) * quotient;
                        println!("      Equation: α^i * (p(z) - p(x)) * (z - x)^(-1)");
                        println!(
                            "      contribution = {:?} * ({:?} - {:?}) * {:?}",
                            alpha_pow, p_at_z, p_at_x, quotient
                        );
                        println!("      contribution = {:?}", contribution);

                        // Compute the reduced opening contribution
                        *ro += contribution;
                        println!("      Updated reduced_opening[{}] = {:?}", log_height, ro);

                        *alpha_pow *= alpha;
                        println!("      Updated α^i = {:?}", alpha_pow);
                    }
                }
            }
        }

        // VALIDATION: Check for constant matrices (lines 429-433 in open_input)
        // If there's a matrix of height 1 (constant), the reduced opening must be zero
        if let Some((_, ro)) = reduced_openings_map.get(&log_blowup) {
            if !ro.is_zero() {
                return Err(
                    "FinalPolyMismatch: Constant matrix has non-zero reduced opening".into(),
                );
            }
            println!("✓ Constant matrix validation passed");
        }

        // Convert to the format expected by verify_query (descending by log_height)
        let reduced_openings: Vec<(usize, Challenge)> = reduced_openings_map
            .into_iter()
            .rev()
            .map(|(log_height, (_, ro))| (log_height, ro))
            .collect();

        println!("Real reduced openings computed from input proof:");
        for (height, evaluation) in &reduced_openings {
            println!(
                "  Height {} (domain size {}): {:?}",
                height,
                1 << height,
                evaluation
            );
        }

        println!("Input opening verification completed");

        // Step 5c: Implement actual verify_query - the FRI folding chain with real values
        println!("\n=== VERIFY_QUERY SIMULATION - ACTUAL FRI FOLDING ===");

        let mut domain_index = index; // Copy for folding simulation

        // VALIDATION: Check reduced openings are valid (lines 244-246 in verify_query)
        let mut ro_iter = reduced_openings.into_iter().peekable();
        if ro_iter.peek().is_none() {
            return Err(
                "InvalidProofShape: No reduced openings - committed to no polynomials".into(),
            );
        }
        if ro_iter.peek().unwrap().0 != log_max_height {
            return Err(format!(
                "InvalidProofShape: First reduced opening height {} != max height {}",
                ro_iter.peek().unwrap().0,
                log_max_height
            )
            .into());
        }
        println!("✓ Reduced openings validation passed");

        // Start with the initial reduced opening (lines 244-247 in verify_query)
        let mut folded_eval = ro_iter.next().unwrap().1;
        println!(
            "Initial folded_eval from reduced openings: {:?}",
            folded_eval
        );

        // Perform the actual FRI folding steps (lines 251-298 in verify_query)
        for (phase_idx, ((beta, commit), opening)) in betas
            .iter()
            .zip(commit_phase_commits.iter())
            .zip(query_proof.commit_phase_openings.iter())
            .enumerate()
        {
            let log_folded_height = log_max_height - phase_idx - 1;

            println!("\n--- Phase {} FRI Folding ---", phase_idx);
            println!("  Beta: {:?}", beta);
            println!("  Commitment: {:?}", commit);
            println!("  Sibling value: {:?}", opening.sibling_value);
            println!(
                "  Domain height: {} (size: {})",
                log_folded_height,
                1 << log_folded_height
            );
            println!("  Current folded_eval: {:?}", folded_eval);

            // Get the index of the other sibling (line 258)
            let index_sibling = domain_index ^ 1;
            println!(
                "  Current index: {}, Sibling index: {}",
                domain_index, index_sibling
            );

            // Create the pair of evaluations (lines 260-261)
            let mut evals = vec![folded_eval; 2];
            evals[index_sibling % 2] = opening.sibling_value;
            println!("  Evaluation pair: {:?}", evals);

            // CRITICAL: REAL MMCS verification of commit phase opening (lines 272-280 in verifier.rs)
            let dims = &[p3_matrix::Dimensions {
                width: 2,
                height: 1 << log_folded_height,
            }];
            println!("  MMCS verification dimensions: {:?}", dims);

            // Update index to parent (line 269)
            domain_index >>= 1;
            println!("  Parent index after fold: {}", domain_index);

            // REAL MMCS verification call using the actual FRI MMCS from PCS
            println!("  === MMCS COMMIT PHASE VERIFICATION IMMEDIATE VALUES ===");
            println!("  Function: fri_mmcs.verify_batch()");
            println!("  Inputs:");
            println!("    commit: {:?}", commit);
            println!("    dims: {:?}", dims);
            println!("    domain_index: {}", domain_index);
            println!("    evaluations: {:?}", evals);
            println!("    opening_proof length: {}", opening.opening_proof.len());
            for (i, proof_element) in opening.opening_proof.iter().enumerate() {
                println!("      proof[{}]: {:?}", i, proof_element);
            }
            println!("  Equation: verify_batch(C_i, dims_i, idx_parent, π_sibling) → bool");
            println!("    where C_i = commit phase i commitment");
            println!("          dims_i = [width=2, height=2^(h-i)]");
            println!("          idx_parent = domain_index >> 1");
            println!("          π_sibling = proof of sibling evaluation");

            // ACTUAL MMCS verification call with the REAL FRI MMCS
            let evals_slice = [evals.clone()];
            let batch_opening_ref =
                p3_commit::BatchOpeningRef::new(&evals_slice, &opening.opening_proof);
            let verification_result =
                pcs.fri_params()
                    .mmcs
                    .verify_batch(commit, dims, domain_index, batch_opening_ref);

            match verification_result {
                Ok(()) => {
                    println!("  Output: SUCCESS (verification passed)");
                    println!("  ✓ MMCS commit phase verification PASSED - sibling proof verified");
                }
                Err(e) => {
                    println!("  Output: FAILED - {:?}", e);
                    return Err(
                        format!("CommitPhaseMmcsError: MMCS verification failed: {:?}", e).into(),
                    );
                }
            }

            // Perform the actual folding operation (line 283)
            println!("  === FRI FOLDING OPERATION IMMEDIATE VALUES ===");
            println!("  Function: TwoAdicFriFolding.fold_row()");
            println!(
                "  This implements the real TwoAdicFriFolding.fold_row logic from two_adic_pcs.rs:106-132"
            );

            let (e0, e1) = (evals[0], evals[1]);
            println!("  Input evaluations:");
            println!("    e0 = evals[0] = {:?}", e0);
            println!("    e1 = evals[1] = {:?}", e1);

            // Compute the subgroup points xs[0] and xs[1] where the evaluations are taken
            let log_arity = 1;
            println!("  Subgroup point computation:");
            println!("    log_arity = {}", log_arity);
            println!("    log_folded_height = {}", log_folded_height);
            println!("    domain_index = {}", domain_index);

            let rev_bits = p3_util::reverse_bits_len(domain_index, log_folded_height);
            println!(
                "    reverse_bits_len({}, {}) = {}",
                domain_index, log_folded_height, rev_bits
            );

            let generator = Val::two_adic_generator(log_folded_height + log_arity);
            println!(
                "    generator = two_adic_generator({}) = {}",
                log_folded_height + log_arity,
                generator
            );

            let subgroup_start = generator.exp_u64(rev_bits as u64);
            println!(
                "    subgroup_start = {}^{} = {}",
                generator, rev_bits, subgroup_start
            );

            let x0 = subgroup_start;
            let x1 = -subgroup_start; // xs[1] = -xs[0] for arity=2
            println!("    x0 = subgroup_start = {}", x0);
            println!("    x1 = -subgroup_start = {}", x1);

            // Real folding: interpolate between (x0, e0) and (x1, e1), evaluate at beta
            println!("  Interpolation formula:");
            println!("    f(X) = e0 + (X - x0) * (e1 - e0) * (x1 - x0)^(-1)");
            println!("    folded_eval = f(β)");

            let beta_minus_x0 = *beta - x0;
            let e1_minus_e0 = e1 - e0;
            let x1_minus_x0 = x1 - x0;
            let inverse = x1_minus_x0.inverse();

            println!("  Step-by-step computation:");
            println!("    β - x0 = {:?} - {} = {:?}", beta, x0, beta_minus_x0);
            println!("    e1 - e0 = {:?} - {:?} = {:?}", e1, e0, e1_minus_e0);
            println!("    x1 - x0 = {} - {} = {}", x1, x0, x1_minus_x0);
            println!("    (x1 - x0)^(-1) = {}^(-1) = {}", x1_minus_x0, inverse);

            let intermediate = beta_minus_x0 * e1_minus_e0 * inverse;
            println!(
                "    (β - x0) * (e1 - e0) * (x1 - x0)^(-1) = {:?}",
                intermediate
            );

            folded_eval = e0 + intermediate;
            println!("  Final result:");
            println!(
                "    folded_eval = e0 + intermediate = {:?} + {:?} = {:?}",
                e0, intermediate, folded_eval
            );

            // Check for roll-ins at this height (lines 295-297)
            println!("  === ROLL-IN CHECK ===");
            if let Some((_, ro)) = ro_iter.next_if(|(lh, _)| *lh == log_folded_height) {
                println!("  Roll-in detected at height {}", log_folded_height);
                println!("  Roll-in equation: folded_eval += β² * reduced_opening");

                let beta_squared = beta.square();
                println!("    β² = ({:?})² = {:?}", beta, beta_squared);
                println!("    reduced_opening = {:?}", ro);

                let roll_in = beta_squared * ro;
                println!(
                    "    roll_in = β² * ro = {:?} * {:?} = {:?}",
                    beta_squared, ro, roll_in
                );

                let old_folded_eval = folded_eval;
                folded_eval += roll_in;
                println!(
                    "    folded_eval = {} + {} = {:?}",
                    old_folded_eval, roll_in, folded_eval
                );
            } else {
                println!("  No roll-in at height {}", log_folded_height);
            }

            println!("  Final folded_eval for this phase: {:?}", folded_eval);
        }

        // VALIDATION: Check if all reduced openings were consumed (lines 301-303 in verify_query)
        if ro_iter.next().is_some() {
            return Err("InvalidProofShape: Failed to fold in some polynomial evaluations".into());
        }
        println!("✓ All reduced openings consumed during folding");

        println!(
            "\nFinal folded evaluation after all FRI folds: {:?}",
            folded_eval
        );

        // Step 5d: Final polynomial evaluation (lines 173-183)
        println!("\n=== FINAL POLYNOMIAL EVALUATION ===");
        println!("Computing final polynomial evaluation at domain point");

        let rev_domain_index = p3_util::reverse_bits_len(domain_index, log_max_height);
        println!(
            "  reverse_bits_len({}, {}) = {}",
            domain_index, log_max_height, rev_domain_index
        );

        let generator = Val::two_adic_generator(log_max_height);
        println!(
            "  generator = two_adic_generator({}) = {}",
            log_max_height, generator
        );

        let x = generator.exp_u64(rev_domain_index as u64);
        println!("  x = {}^{} = {}", generator, rev_domain_index, x);
        println!("  Final evaluation point x: {}", x);
        println!("  Final domain index: {}", domain_index);

        // Evaluate final polynomial - simplified for constant polynomial (log_final_poly_len=0)
        println!("  Final polynomial evaluation (simplified for constant):");
        println!("    Final polynomial is constant: {:?}", final_poly[0]);
        println!("    No Horner's method needed - polynomial evaluation = constant");
        let eval = final_poly[0]; // Constant polynomial evaluation
        println!("  Final polynomial evaluation: {:?}", eval);

        // CRITICAL: Final polynomial mismatch check (matching verifier.rs:185-187)
        println!("\n=== FINAL POLYNOMIAL MISMATCH CHECK ===");
        println!("Final polynomial evaluation at x: {:?}", eval);
        println!("Folded evaluation from FRI process: {:?}", folded_eval);

        if eval != folded_eval {
            println!("ERROR: Final polynomial mismatch detected!");
            println!("  Expected (final poly): {:?}", eval);
            println!("  Actual (folded eval):  {:?}", folded_eval);
            return Err("FinalPolyMismatch".into());
        }

        println!("✓ Final polynomial evaluation matches folded evaluation");
        println!("  Both evaluate to: {:?}", eval);

        println!("Final verification check simulation completed");
    }

    Ok(())
}
