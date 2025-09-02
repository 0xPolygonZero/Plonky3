use itertools::Itertools;
use itertools::zip_eq;
use p3_challenger::CanObserve;
use p3_challenger::CanSampleBits;
use p3_challenger::FieldChallenger;
use p3_commit::Pcs;
use p3_commit::PolynomialSpace;
use p3_field::BasedVectorSpace;
use p3_field::ExtensionField;
use p3_field::Field;
use p3_field::PrimeCharacteristicRing;
use p3_field::extension::BinomiallyExtendable;
use p3_uni_stark::Domain;
use p3_uni_stark::Proof;
use p3_uni_stark::StarkGenericConfig;
use p3_uni_stark::Val;

use crate::circuit_builder::CircuitError;
use crate::circuit_builder::ExtensionWireId;
use crate::circuit_builder::gates::arith_gates::AddExtensionGate;
use crate::circuit_builder::gates::arith_gates::MulExtensionGate;
use crate::circuit_builder::gates::arith_gates::SubExtensionGate;
use crate::circuit_builder::{CircuitBuilder, WireId};
use crate::verifier::recursive_traits::CommitForRecursiveVerif;
use crate::verifier::recursive_traits::CommitRecursiveVerif;
use crate::verifier::recursive_traits::PcsRecursiveVerif;
use crate::verifier::recursive_traits::RecursiveAir;
use crate::verifier::recursive_traits::RecursiveStarkGenerationConfig;

#[derive(Clone)]
pub struct FriProofWires<const D: usize, Comm: CommitRecursiveVerif, InputProof> {
    pub commit_phase_commits: Vec<Comm>,
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
pub struct CommitmentWires<Comm: CommitRecursiveVerif> {
    pub trace_wires: Comm,
    pub quotient_chunks_wires: Comm,
    pub random_commit: Option<Comm>,
}

// TODO: Move these structures to their respective crates.
#[derive(Clone)]
pub struct OpenedValuesWires<const D: usize> {
    pub trace_local_wires: Vec<ExtensionWireId<D>>,
    pub trace_next_wires: Vec<ExtensionWireId<D>>,
    pub quotient_chunks_wires: Vec<Vec<ExtensionWireId<D>>>,
    pub random: Option<Vec<ExtensionWireId<D>>>,
}

// Simulated method as we don't have the sponges yet.
pub fn generate_challenges<
    SC: StarkGenericConfig,
    Comm: CommitRecursiveVerif,
    InputProof,
    const D: usize,
    const DIGEST_ELEMS: usize,
>(
    config: &SC,
    circuit: &mut CircuitBuilder<Val<SC>, D>,
    proof_wires: &ProofWires<D, Comm, InputProof>,
    proof: &Proof<SC>,
    public_values: &Vec<WireId>,
    challenges: &Vec<ExtensionWireId<D>>,
    init_trace_domain: &Domain<SC>,
) where
    <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment: CommitForRecursiveVerif<Val<SC>>,
{
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
        degree_bits: _,
    } = proof_wires;

    let Proof {
        commitments,
        opened_values: _,
        opening_proof: _,
        degree_bits,
    } = proof;

    challenger.observe(Val::<SC>::from_usize(*degree_bits));
    challenger.observe(Val::<SC>::from_usize(*degree_bits - config.is_zk()));

    let local_comm = commitments.trace.get_values();
    let local_comm_wires = <Comm as CommitRecursiveVerif>::get_wires(trace_wires);

    // for v in local_wire_vals {
    //     challenger.observe(v.as_base().expect("The output should be a basis element"));
    // }
    challenger.observe_slice(&local_comm);
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

    let quotient_chunks_comm_wires = quotient_chunks_wires.get_wires();
    let quotient_chunks_vals = commitments.quotient_chunks.get_values();
    circuit
        .set_wire_values(&quotient_chunks_comm_wires, quotient_chunks_vals.as_slice())
        .unwrap();

    // let quotient_chunk_com_vals = quotient_chunks_wires
    //     .iter()
    //     .map(|wire| {
    //         circuit
    //             .get_wire_value(*wire)
    //             .unwrap()
    //             .unwrap()
    //             .as_base()
    //             .unwrap()
    //     })
    //     .collect::<Vec<_>>();
    challenger.observe_slice(&quotient_chunks_vals);

    if let Some(r_commit) = commitments.random.as_ref() {
        let r_commit_vals = r_commit.get_values();
        // let r_commit_vals = r_commit
        //     .iter()
        //     .map(|wire| {
        //         circuit
        //             .get_wire_value(*wire)
        //             .unwrap()
        //             .unwrap()
        //             .as_base()
        //             .unwrap()
        //     })
        //     .collect::<Vec<_>>();
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
            .get_wires()
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
        // let comm = comm_wires
        //     .iter()
        //     .map(|wire| {
        //         circuit
        //             .get_wire_value(*wire)
        //             .unwrap()
        //             .unwrap()
        //             .as_base()
        //             .unwrap()
        //     })
        //     .collect::<Vec<_>>();
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
        circuit.set_extension_wires(*c, challenge_vals).unwrap();
    }
}

fn get_circuit_challenges<
    SC: RecursiveStarkGenerationConfig<InputProof, D>,
    const D: usize,
    InputProof,
>(
    proof_wires: &ProofWires<D, SC::Comm, InputProof>,
    circuit: &mut CircuitBuilder<SC::Val, D>,
) -> Vec<ExtensionWireId<D>> {
    let mut challenges = vec![];
    // Observe degree bits and degree_bits - is_zk.
    // Observe local wires.
    // Observe public values.
    challenges.push(circuit.new_extension_wires());
    // Observe quotient chunks.
    // Observe random commitment if any.
    // zeta and zeta_next
    challenges.push(circuit.new_extension_wires());
    challenges.push(circuit.new_extension_wires());

    let pcs_challenges =
        <SC as RecursiveStarkGenerationConfig<InputProof, D>>::Pcs::get_challenges_circuit(
            circuit,
            proof_wires,
        );

    challenges.extend(pcs_challenges);

    challenges
}

#[derive(Clone)]
pub struct ProofWires<const D: usize, Comm: CommitRecursiveVerif, InputProof> {
    pub commitments_wires: CommitmentWires<Comm>,
    pub opened_values_wires: OpenedValuesWires<D>,
    pub fri_proof: FriProofWires<D, Comm, InputProof>,
    degree_bits: usize,
}

pub fn verify_circuit<
    A,
    SC: RecursiveStarkGenerationConfig<InputProof, D>,
    InputProof,
    const D: usize,
    const DIGEST_ELEMS: usize,
>(
    config: &SC,
    air: &A,
    proof_wires: &ProofWires<D, SC::Comm, InputProof>,
    public_values: &Vec<WireId>,
) -> Result<(), CircuitError>
where
    SC::Val: BinomiallyExtendable<D>,
    A: RecursiveAir<SC::Val, D>,
    InputProof: Clone,
    <SC as RecursiveStarkGenerationConfig<InputProof, D>>::Pcs:
        PcsRecursiveVerif<InputProof, SC::Comm, SC::Domain, SC::Val, SC::Challenge, D>,
{
    let ProofWires {
        commitments_wires:
            CommitmentWires {
                trace_wires,
                quotient_chunks_wires,
                random_commit,
            },
        opened_values_wires:
            OpenedValuesWires {
                trace_local_wires: opened_trace_local_wires,
                trace_next_wires: opened_trace_next_wires,
                quotient_chunks_wires: opened_quotient_chunks_wires,
                random: opened_random,
            },
        fri_proof:
            FriProofWires {
                commit_phase_commits: _,
                query_proofs: _,
                final_poly: _,
                pow_witness: _,
            },
        degree_bits,
    } = proof_wires;
    let degree = 1 << degree_bits;
    let log_quotient_degree =
        A::get_log_quotient_degree(air, 0, public_values.len(), config.is_zk());
    let quotient_degree = 1 << (log_quotient_degree + config.is_zk());

    let pcs = config.pcs();
    let trace_domain = pcs.natural_domain_for_degree(degree);
    let init_trace_domain = pcs.natural_domain_for_degree(degree >> (config.is_zk()));
    let mut circuit = CircuitBuilder::<SC::Val, D>::new();

    let quotient_domain =
        pcs.create_disjoint_domain(trace_domain, 1 << (degree_bits + log_quotient_degree));
    let quotient_chunks_domains = pcs.split_domains(&quotient_domain, quotient_degree);

    let randomized_quotient_chunks_domains = quotient_chunks_domains
        .iter()
        .map(|domain| pcs.natural_domain_for_degree(pcs.size(domain) << (config.is_zk())))
        .collect_vec();

    // Challenger is called here. But we don't have the interactions or hash tables yet.
    // TODO: We need to simulate it for now.
    let challenge_wires = get_circuit_challenges::<SC, D, InputProof>(proof_wires, &mut circuit);

    // Verify shape.
    let air_width = A::width(air);
    let validate_shape = opened_trace_local_wires.len() == air_width
        && opened_trace_next_wires.len() == air_width
        && opened_quotient_chunks_wires.len() == quotient_degree
        && opened_quotient_chunks_wires
            .iter()
            .all(|opened_chunk| opened_chunk.len() == SC::Challenge::DIMENSION);
    if !validate_shape {
        return Err(CircuitError::InvalidProofShape);
    }

    let alpha: ExtensionWireId<D> = challenge_wires[0];
    let zeta: ExtensionWireId<D> = challenge_wires[1];
    let zeta_next: ExtensionWireId<D> = challenge_wires[2];

    // Need to simulate Fri here.
    let mut coms_to_verify = if let Some(r_commit) = &random_commit {
        let random_values = opened_random
            .as_ref()
            .ok_or(CircuitError::RandomizationError)?;
        vec![(
            r_commit,
            vec![(trace_domain, vec![(zeta, random_values.clone())])],
        )]
    } else {
        vec![]
    };
    coms_to_verify.extend(vec![
        (
            trace_wires,
            vec![(
                trace_domain,
                vec![
                    (zeta, opened_trace_local_wires.clone()),
                    (zeta_next, opened_trace_next_wires.clone()),
                ],
            )],
        ),
        (
            quotient_chunks_wires,
            // Check the commitment on the randomized domains.
            zip_eq(
                randomized_quotient_chunks_domains.iter(),
                opened_quotient_chunks_wires,
            )
            .map(|(domain, values)| (*domain, vec![(zeta, values.clone())]))
            .collect_vec(),
        ),
    ]);
    pcs.verify_circuit(&mut circuit, &challenge_wires[3..], &coms_to_verify);

    let zero = circuit.add_extension_constant(SC::Challenge::ZERO);
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
                    let v_n = vanishing_poly_at_point_circuit::<
                        InputProof,
                        SC::Comm,
                        SC,
                        SC::Domain,
                        D,
                    >(config, *other_domain, zeta, &mut circuit);

                    let first_point = circuit
                        .add_extension_constant(SC::Challenge::from(pcs.first_point(domain)));
                    let other_v_n =
                        vanishing_poly_at_point_circuit::<InputProof, SC::Comm, SC, SC::Domain, D>(
                            config,
                            *other_domain,
                            first_point,
                            &mut circuit,
                        );
                    let div = circuit.new_extension_wires();
                    MulExtensionGate::<SC::Val, D>::add_to_circuit(
                        &mut circuit,
                        other_v_n,
                        div,
                        v_n,
                    );

                    let new_total = circuit.new_extension_wires();
                    MulExtensionGate::<SC::Val, D>::add_to_circuit(
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
                circuit.add_extension_constant(SC::Challenge::ith_basis_element(e_i).unwrap());
            let inner_mul = circuit.new_extension_wires();
            MulExtensionGate::<SC::Val, D>::add_to_circuit(&mut circuit, e_i_wire, *c, inner_mul);
            let new_s = circuit.new_extension_wires();
            AddExtensionGate::<SC::Val, D>::add_to_circuit(&mut circuit, cur_s, inner_mul, new_s);
            cur_s = inner_mul;
        }
        let mul = circuit.new_extension_wires();
        MulExtensionGate::<SC::Val, D>::add_to_circuit(&mut circuit, cur_s, zp, mul);
        let add_wire = circuit.new_extension_wires();
        AddExtensionGate::<SC::Val, D>::add_to_circuit(&mut circuit, quotient, mul, add_wire);
        quotient = add_wire;
    }

    let sels = pcs.selectors_at_point_circuit(&mut circuit, &init_trace_domain, &zeta);
    let folded_constraints = air.eval_folded_circuit::<SC::Challenge>(
        &mut circuit,
        &sels,
        &alpha,
        &vec![],
        &vec![],
        opened_trace_local_wires,
        opened_trace_next_wires,
        &public_values,
    );

    // Compute folded_constraints * sels.inv_vanishing.
    let folded_mul = circuit.new_extension_wires();
    MulExtensionGate::<SC::Val, D>::add_to_circuit(
        &mut circuit,
        folded_constraints,
        sels.inv_vanishing,
        folded_mul,
    );

    // Check that folded_constraints * sels.inv_vanishing == quotient
    SubExtensionGate::<SC::Val, D>::add_to_circuit(&mut circuit, folded_mul, quotient, zero);

    Ok(())
}

fn vanishing_poly_at_point_circuit<
    InputProof,
    Comm: CommitRecursiveVerif,
    SC: RecursiveStarkGenerationConfig<InputProof, D>,
    Domain,
    const D: usize,
>(
    config: &SC,
    domain: Domain,
    zeta: ExtensionWireId<D>,
    circuit: &mut CircuitBuilder<SC::Val, D>,
) -> ExtensionWireId<D>
where
    SC::Val: BinomiallyExtendable<D>,
    <SC as RecursiveStarkGenerationConfig<InputProof, D>>::Pcs:
        PcsRecursiveVerif<InputProof, Comm, Domain, SC::Val, SC::Challenge, D>,
{
    let pcs = config.pcs();
    let inv =
        circuit.add_extension_constant(SC::Challenge::from(pcs.first_point(&domain).inverse()));

    let mul = circuit.new_extension_wires();
    MulExtensionGate::<SC::Val, D>::add_to_circuit(circuit, zeta, inv, mul);
    let size_wire = circuit.add_extension_constant(SC::Challenge::from_usize(pcs.size(&domain)));
    let exp = circuit.new_extension_wires();
    MulExtensionGate::<SC::Val, D>::add_to_circuit(circuit, mul, size_wire, exp);

    let one = circuit.add_extension_constant(SC::Challenge::ONE);
    let v_n = circuit.new_extension_wires();
    SubExtensionGate::<SC::Val, D>::add_to_circuit(circuit, exp, one, v_n);

    v_n
}
