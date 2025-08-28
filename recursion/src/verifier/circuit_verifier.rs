use itertools::Itertools;
use p3_air::Air;
use p3_challenger::CanObserve;
use p3_challenger::CanSampleBits;
use p3_challenger::FieldChallenger;
use p3_commit::{Pcs, PolynomialSpace};
use p3_field::BasedVectorSpace;
use p3_field::ExtensionField;
use p3_field::Field;
use p3_field::PrimeCharacteristicRing;
use p3_field::extension::BinomiallyExtendable;

use p3_uni_stark::Domain;
use p3_uni_stark::get_log_quotient_degree;
use p3_uni_stark::{PcsError, StarkGenericConfig, SymbolicAirBuilder, Val, VerificationError};

use crate::circuit_builder::ChallengeWireId;
use crate::circuit_builder::gates::arith_gates::AddExtensionGate;
use crate::circuit_builder::gates::arith_gates::MulExtensionGate;
use crate::circuit_builder::gates::arith_gates::SubExtensionGate;
use crate::circuit_builder::{CircuitBuilder, WireId};

#[derive(Clone)]
pub struct FriProofWires<const DIGEST_ELEMS: usize, InputProof> {
    pub commit_phase_commits: Vec<[usize; DIGEST_ELEMS]>,
    pub query_proofs: Vec<QueryProofWires<InputProof>>,
    pub final_poly: Vec<usize>,
    pub pow_witness: usize,
}

#[derive(Clone)]
pub struct QueryProofWires<InputProof> {
    pub input_proof: InputProof,
    pub commit_phase_openings: Vec<CommitPhaseProofStepWires>,
}

#[derive(Clone)]
pub struct CommitPhaseProofStepWires {
    pub sibling_value: usize,
    pub opening_proof: Vec<usize>,
}

#[derive(Clone)]
pub struct CommitmentWires<const DIGEST_ELEMS: usize> {
    pub trace_wires: [WireId; DIGEST_ELEMS],
    pub quotient_chunks_wires: [WireId; DIGEST_ELEMS],
    pub random_commit: Option<[WireId; DIGEST_ELEMS]>,
}

// TODO: Move these structures to their respective crates.
#[derive(Clone)]
pub struct OpenedValuesWires<const D: usize> {
    pub trace_local_wires: Vec<ChallengeWireId<D>>,
    pub trace_next_wires: Vec<ChallengeWireId<D>>,
    pub quotient_chunks_wires: Vec<Vec<ChallengeWireId<D>>>,
    pub random: Option<Vec<ChallengeWireId<D>>>,
}

// Simulated method as we don't have the sponges yet.
pub fn generate_challenges<
    SC: StarkGenericConfig,
    InputProof,
    const D: usize,
    const DIGEST_ELEMS: usize,
>(
    config: &SC,
    circuit: &mut CircuitBuilder<Val<SC>, D>,
    proof_wires: &ProofWires<D, DIGEST_ELEMS, InputProof>,
    public_values: &Vec<WireId>,
    challenges: &Vec<ChallengeWireId<D>>,
    init_trace_domain: &Domain<SC>,
) {
    // This assumes that the public values and proof inputs are already set up in the circuit.
    let mut challenger = config.initialise_challenger();
    let mut challenge_values: Vec<SC::Challenge> = vec![];
    let ProofWires {
        commitments_wires:
            CommitmentWires {
                trace_wires,
                quotient_chunks_wires,
                random_commit,
            },
        opened_values_wires:
            OpenedValuesWires {
                trace_local_wires,
                trace_next_wires,
                quotient_chunks_wires: opened_quotient_chunks_wires,
                random,
            },
        fri_proof:
            FriProofWires {
                commit_phase_commits,
                query_proofs,
                final_poly,
                pow_witness: _,
            },
        degree_bits,
    } = proof_wires;

    challenger.observe(Val::<SC>::from_usize(*degree_bits));
    challenger.observe(Val::<SC>::from_usize(*degree_bits - config.is_zk()));

    let local_wire_vals = trace_wires
        .iter()
        .map(|wire| circuit.get_wire_value(*wire).unwrap().unwrap());
    for v in local_wire_vals {
        challenger.observe(v.as_base().expect("The output should be a basis element"));
    }

    let public_values_vals = public_values
        .iter()
        .map(|wire| {
            circuit
                .get_wire_value(*wire)
                .unwrap()
                .unwrap()
                .as_base()
                .unwrap()
        })
        .collect::<Vec<_>>();
    challenger.observe_slice(&public_values_vals);
    // 1
    challenge_values.push(challenger.sample_algebra_element());

    let quotient_chunk_com_vals = quotient_chunks_wires
        .iter()
        .map(|wire| {
            circuit
                .get_wire_value(*wire)
                .unwrap()
                .unwrap()
                .as_base()
                .unwrap()
        })
        .collect::<Vec<_>>();
    challenger.observe_slice(&quotient_chunk_com_vals);

    if let Some(r_commit) = random_commit {
        let r_commit_vals = r_commit
            .iter()
            .map(|wire| {
                circuit
                    .get_wire_value(*wire)
                    .unwrap()
                    .unwrap()
                    .as_base()
                    .unwrap()
            })
            .collect::<Vec<_>>();
        challenger.observe_slice(&r_commit_vals);
    }

    let zeta = challenger.sample_algebra_element();
    // 2 and 3
    challenge_values.push(zeta);
    let zeta_next = init_trace_domain.next_point(zeta).unwrap();
    challenge_values.push(zeta_next);

    let mut coms_to_verify = if let Some(r_commit) = &random_commit {
        let random_values = random.as_ref().unwrap();
        vec![(r_commit.clone(), vec![vec![(zeta, random_values.clone())]])]
    } else {
        vec![]
    };
    coms_to_verify.extend(vec![
        (
            trace_wires.clone(),
            vec![vec![
                (zeta, trace_local_wires.clone()),
                (zeta_next, trace_next_wires.clone()),
            ]],
        ),
        (
            quotient_chunks_wires.clone(),
            // Check the commitment on the randomized domains.
            opened_quotient_chunks_wires
                .iter()
                .map(|w| vec![(zeta, w.clone())])
                .collect::<Vec<_>>(),
        ),
    ]);

    for (_, round) in &coms_to_verify {
        for mat in round {
            for (_, point) in mat {
                point.iter().for_each(|&opening| {
                    for o in opening {
                        challenger.observe(
                            circuit
                                .get_wire_value(o)
                                .unwrap()
                                .unwrap()
                                .as_base()
                                .unwrap(),
                        )
                    }
                })
            }
        }
    }

    // FRI challenges.
    // Batch combination challenge.
    challenge_values.push(challenger.sample_algebra_element());

    // betas
    commit_phase_commits.iter().for_each(|comm_wires| {
        let comm = comm_wires
            .iter()
            .map(|wire| {
                circuit
                    .get_wire_value(*wire)
                    .unwrap()
                    .unwrap()
                    .as_base()
                    .unwrap()
            })
            .collect::<Vec<_>>();
        challenger.observe_slice(&comm);
        challenge_values.push(challenger.sample_algebra_element());
    });

    final_poly.iter().for_each(|x| {
        let value = circuit
            .get_wire_value(*x)
            .unwrap()
            .unwrap()
            .as_base()
            .unwrap();
        challenger.observe(value);
    });

    let (log_blowup, log_final_poly_len) = config.pcs().get_log_blowup_final_height();
    let log_max_height = commit_phase_commits.len() + log_blowup + log_final_poly_len;

    for _query_proof in query_proofs {
        challenge_values.push(SC::Challenge::from_usize(
            // extra_query_bits = 0 for two-adic pcs.
            challenger.sample_bits(log_max_height),
        ));
    }

    for (c, c_v) in challenges.iter().zip_eq(challenge_values.iter()) {
        let challenge_vals = c_v.as_basis_coefficients_slice().try_into().unwrap();
        circuit.set_challenge_wires(*c, challenge_vals).unwrap();
    }
}

fn get_circuit_challenges<
    SC: StarkGenericConfig,
    const D: usize,
    const DIGEST_ELEMS: usize,
    InputProof,
>(
    proof_wires: &ProofWires<D, DIGEST_ELEMS, InputProof>,
    circuit: &mut CircuitBuilder<Val<SC>, D>,
) -> Vec<ChallengeWireId<D>> {
    // Observe degree bits and degree_bits - is_zk.
    // Observe local wires.
    // Observe public values.
    let mut num_challenges = 1; // alpha
    // Observe quotient chunks.
    // Observe random commitment if any.
    num_challenges += 2; // zeta and zeta_next
    // Observe coms_to_verify.
    num_challenges += 1; // batch combination challenge

    // Observe the betas at the same time
    num_challenges += proof_wires.fri_proof.commit_phase_commits.len(); // betas

    // Observe final poly.
    num_challenges += proof_wires.fri_proof.query_proofs.len(); // final poly evaluations

    let mut challenges = Vec::with_capacity(num_challenges);
    for _ in 0..num_challenges {
        challenges.push(circuit.new_challenge_wires());
    }

    challenges
}

#[derive(Clone)]
pub struct ProofWires<const D: usize, const DIGEST_ELEMS: usize, InputProof> {
    pub commitments_wires: CommitmentWires<DIGEST_ELEMS>,
    pub opened_values_wires: OpenedValuesWires<D>,
    pub fri_proof: FriProofWires<DIGEST_ELEMS, InputProof>,
    degree_bits: usize,
}

pub fn verify_circuit<
    A,
    SC: StarkGenericConfig,
    InputProof,
    const D: usize,
    const DIGEST_ELEMS: usize,
>(
    config: &SC,
    air: &A,
    proof_wires: ProofWires<D, DIGEST_ELEMS, InputProof>,
    public_values: &Vec<WireId>,
) -> Result<(), VerificationError<PcsError<SC>>>
where
    Val<SC>: BinomiallyExtendable<D>,
    A: Air<SymbolicAirBuilder<Val<SC>>>,
    InputProof: Clone,
{
    let ProofWires {
        commitments_wires:
            CommitmentWires {
                trace_wires: _,
                quotient_chunks_wires,
                random_commit: _,
            },
        opened_values_wires:
            OpenedValuesWires {
                trace_local_wires,
                trace_next_wires,
                quotient_chunks_wires: opened_quotient_chunks_wires,
                random: _,
            },
        fri_proof:
            FriProofWires {
                commit_phase_commits: _,
                query_proofs: _,
                final_poly: _,
                pow_witness: _,
            },
        degree_bits,
    } = proof_wires.clone();
    let degree = 1 << degree_bits;
    let log_quotient_degree =
        get_log_quotient_degree::<Val<SC>, A>(air, 0, public_values.len(), config.is_zk());
    let quotient_degree = 1 << (log_quotient_degree + config.is_zk());

    let pcs = config.pcs();
    let trace_domain = pcs.natural_domain_for_degree(degree);
    let init_trace_domain = pcs.natural_domain_for_degree(degree >> (config.is_zk()));
    let mut circuit = CircuitBuilder::<Val<SC>, D>::new();

    let quotient_domain =
        trace_domain.create_disjoint_domain(1 << (degree_bits + log_quotient_degree));
    let quotient_chunks_domains = quotient_domain.split_domains(quotient_degree);

    // Challenger is called here. But we don't have the interactions or hash tables yet.
    // So I need to simulate it.
    let challenge_wires =
        get_circuit_challenges::<SC, D, DIGEST_ELEMS, InputProof>(&proof_wires, &mut circuit);
    generate_challenges(
        config,
        &mut circuit,
        &proof_wires,
        public_values,
        &challenge_wires,
        &init_trace_domain,
    );

    let air_width = A::width(air);
    let validate_shape = trace_local_wires.len() == air_width
        && trace_next_wires.len() == air_width
        && quotient_chunks_wires.len() == quotient_chunks_wires.len()
        && opened_quotient_chunks_wires
            .iter()
            .all(|opened_chunk| opened_chunk.len() == SC::Challenge::DIMENSION);
    if !validate_shape {
        return Err(VerificationError::InvalidProofShape);
    }

    let alpha: ChallengeWireId<D> = challenge_wires[0];
    let zeta: ChallengeWireId<D> = challenge_wires[1];
    let zeta_next: ChallengeWireId<D> = challenge_wires[2];

    // Need to simulate Fri here.

    let zero = circuit.add_challenge_constant(SC::Challenge::ZERO);
    let zps = quotient_chunks_domains
        .iter()
        .enumerate()
        .map(|(i, domain)| {
            let mut total = zero;
            quotient_chunks_domains
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .for_each(|(_, other_domain)| {
                    let v_n =
                        vanishing_poly_at_point_circuit::<SC, D>(*other_domain, zeta, &mut circuit);

                    let first_point =
                        circuit.add_challenge_constant(SC::Challenge::from(domain.first_point()));
                    let other_v_n = vanishing_poly_at_point_circuit::<SC, D>(
                        *other_domain,
                        first_point,
                        &mut circuit,
                    );
                    let div = circuit.new_challenge_wires();
                    MulExtensionGate::<Val<SC>, D>::add_to_circuit(
                        &mut circuit,
                        other_v_n,
                        div,
                        v_n,
                    );

                    let new_total = circuit.new_challenge_wires();
                    MulExtensionGate::<Val<SC>, D>::add_to_circuit(
                        &mut circuit,
                        total,
                        v_n,
                        new_total,
                    );
                    total = new_total;
                });
            total
        })
        .collect_vec();

    let mut quotient = zero;
    for (i, chunk) in opened_quotient_chunks_wires.iter().enumerate() {
        let zp = zps[i];

        let mut cur_s = zero;
        for (e_i, c) in chunk.iter().enumerate() {
            let e_i_wire =
                circuit.add_challenge_constant(SC::Challenge::ith_basis_element(e_i).unwrap());
            let inner_mul = circuit.new_challenge_wires();
            MulExtensionGate::<Val<SC>, D>::add_to_circuit(&mut circuit, e_i_wire, *c, inner_mul);
            let new_s = circuit.new_challenge_wires();
            AddExtensionGate::<Val<SC>, D>::add_to_circuit(&mut circuit, cur_s, inner_mul, new_s);
            cur_s = inner_mul;
        }
        let mul = circuit.new_challenge_wires();
        MulExtensionGate::<Val<SC>, D>::add_to_circuit(&mut circuit, cur_s, zp, mul);
        let add_wire = circuit.new_challenge_wires();
        AddExtensionGate::<Val<SC>, D>::add_to_circuit(&mut circuit, quotient, mul, add_wire);
        quotient = add_wire;
    }

    // let sels = selectors_at_point_circuit(init_trace_domain, zeta);

    Ok(())
}

fn vanishing_poly_at_point_circuit<SC: StarkGenericConfig, const D: usize>(
    domain: Domain<SC>,
    zeta: ChallengeWireId<D>,
    circuit: &mut CircuitBuilder<Val<SC>, D>,
) -> ChallengeWireId<D>
where
    Val<SC>: BinomiallyExtendable<D>,
{
    let inv = circuit.add_challenge_constant(SC::Challenge::from(domain.first_point().inverse()));

    let mul = circuit.new_challenge_wires();
    MulExtensionGate::<Val<SC>, D>::add_to_circuit(circuit, zeta, inv, mul);
    let size_wire = circuit.add_challenge_constant(SC::Challenge::from_usize(domain.size()));
    let exp = circuit.new_challenge_wires();
    MulExtensionGate::<Val<SC>, D>::add_to_circuit(circuit, mul, size_wire, exp);

    let one = circuit.add_challenge_constant(SC::Challenge::ONE);
    let v_n = circuit.new_challenge_wires();
    SubExtensionGate::<Val<SC>, D>::add_to_circuit(circuit, exp, one, v_n);

    v_n
}
