use itertools::Itertools;
use itertools::zip_eq;
use p3_challenger::CanObserve;
use p3_challenger::FieldChallenger;
use p3_commit::Pcs;
use p3_commit::PolynomialSpace;
use p3_field::BasedVectorSpace;
use p3_field::ExtensionField;
use p3_field::Field;
use p3_field::PrimeCharacteristicRing;
use p3_field::extension::BinomiallyExtendable;
use p3_uni_stark::Commitments;
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
use crate::verifier::recursive_traits::ForRecursiveVersion;
use crate::verifier::recursive_traits::PcsGeneration;
use crate::verifier::recursive_traits::PcsRecursiveVerif;
use crate::verifier::recursive_traits::RecursiveAir;
use crate::verifier::recursive_traits::RecursiveStarkGenerationConfig;
use crate::verifier::recursive_traits::RecursiveVersion;

#[derive(Clone)]
pub struct CommitmentWires<Comm: RecursiveVersion> {
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

// Simulated method to generate actual challenge values as we don't have the sponges yet.
pub fn generate_challenges<
    SC: StarkGenericConfig,
    InputProof: ForRecursiveVersion<Val<SC>>,
    const D: usize,
>(
    config: &SC,
    circuit: &mut CircuitBuilder<Val<SC>, D>,
    proof: &Proof<SC>,
    public_values: &Vec<WireId>,
    challenges: &Vec<ExtensionWireId<D>>,
    init_trace_domain: &Domain<SC>,
) where
    <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment: ForRecursiveVersion<Val<SC>>,
    <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Proof: ForRecursiveVersion<Val<SC>>,
    SC::Pcs:
        PcsRecursiveVerif<
                InputProof::RecVerif,
                <<SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Proof as ForRecursiveVersion<
                    Val<SC>,
                >>::RecVerif,
                <<SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment as ForRecursiveVersion<
                    Val<SC>,
                >>::RecVerif,
                Domain<SC>,
                Val<SC>,
                SC::Challenge,
                D,
            > + PcsGeneration<SC, <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Proof>,
{
    // This assumes that the public values and proof inputs are already set up in the circuit.
    let mut challenger = config.initialise_challenger();
    let mut challenge_values: Vec<SC::Challenge> = vec![];

    let Proof {
        commitments:
            Commitments {
                trace,
                quotient_chunks,
                random,
            },
        opened_values,
        opening_proof,
        degree_bits,
    } = proof;

    challenger.observe(Val::<SC>::from_usize(*degree_bits));
    challenger.observe(Val::<SC>::from_usize(*degree_bits - config.is_zk()));

    let local_comm = trace.get_values();

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
    challenge_values.push(challenger.sample_algebra_element());

    let quotient_chunks_vals = quotient_chunks.get_values();
    challenger.observe_slice(&quotient_chunks_vals);

    if let Some(r_commit) = random.as_ref() {
        let r_commit_vals = r_commit.get_values();
        challenger.observe_slice(&r_commit_vals);
    }

    let zeta = challenger.sample_algebra_element();
    challenge_values.push(zeta);
    let zeta_next = init_trace_domain.next_point(zeta).unwrap();
    challenge_values.push(zeta_next);

    let mut coms_to_verify = if let Some(r_commit) = &random {
        let random_values = opened_values.random.as_ref().unwrap();
        vec![(r_commit.clone(), vec![vec![(zeta, random_values.clone())]])]
    } else {
        vec![]
    };
    coms_to_verify.extend(vec![
        (
            trace.clone(),
            vec![vec![
                (zeta, opened_values.trace_local.clone()),
                (zeta_next, opened_values.trace_next.clone()),
            ]],
        ),
        (
            quotient_chunks.clone(),
            // Check the commitment on the randomized domains.
            opened_values
                .quotient_chunks
                .iter()
                .map(|w| vec![(zeta, w.clone())])
                .collect::<Vec<_>>(),
        ),
    ]);

    for (_, round) in &coms_to_verify {
        for mat in round {
            for (_, point) in mat {
                point.iter().for_each(|&opening| {
                    challenger.observe_algebra_element(opening);
                })
            }
        }
    }

    let pcs_challenges = <SC::Pcs as PcsGeneration<
        SC,
        <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Proof,
    >>::generate_challenges::<InputProof, D>(
        config,
        &mut challenger,
        &coms_to_verify,
        opening_proof,
    );

    challenge_values.extend(pcs_challenges);

    for (c, c_v) in challenges.iter().zip_eq(challenge_values.iter()) {
        let challenge_vals = c_v.as_basis_coefficients_slice().try_into().unwrap();
        circuit.set_extension_wires(*c, challenge_vals).unwrap();
    }
}

// Method to get all the challenge wires.
fn get_circuit_challenges<
    SC: RecursiveStarkGenerationConfig<InputProof, OpeningProof, D>,
    const D: usize,
    InputProof: RecursiveVersion,
    OpeningProof: RecursiveVersion,
>(
    proof_wires: &ProofWires<D, SC::Comm, OpeningProof>,
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
        <SC as RecursiveStarkGenerationConfig<InputProof, OpeningProof, D>>::Pcs::get_challenges_circuit(
            circuit,
            proof_wires,
        );

    challenges.extend(pcs_challenges);

    challenges
}

/// Structure representing all the wires necessary for an input proof.
#[derive(Clone)]
pub struct ProofWires<const D: usize, Comm: RecursiveVersion, OpeningProof: RecursiveVersion> {
    pub commitments_wires: CommitmentWires<Comm>,
    pub opened_values_wires: OpenedValuesWires<D>,
    pub opening_proof: OpeningProof,
    degree_bits: usize,
}

pub fn verify_circuit<
    A,
    SC: RecursiveStarkGenerationConfig<InputProof, OpeningProof, D>,
    InputProof: RecursiveVersion,
    OpeningProof: RecursiveVersion,
    const D: usize,
    const DIGEST_ELEMS: usize,
>(
    config: &SC,
    air: &A,
    proof_wires: &ProofWires<D, SC::Comm, OpeningProof>,
    public_values: &Vec<WireId>,
) -> Result<(), CircuitError>
where
    SC::Val: BinomiallyExtendable<D>,
    A: RecursiveAir<SC::Val, D>,
    InputProof: Clone,
    <SC as RecursiveStarkGenerationConfig<InputProof, OpeningProof, D>>::Pcs: PcsRecursiveVerif<InputProof, OpeningProof, SC::Comm, SC::Domain, SC::Val, SC::Challenge, D>,
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
        opening_proof,
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
    let challenge_wires =
        get_circuit_challenges::<SC, D, InputProof, OpeningProof>(proof_wires, &mut circuit);

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
    pcs.verify_circuit(
        &mut circuit,
        &challenge_wires[3..],
        &coms_to_verify,
        opening_proof,
    );

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
                        OpeningProof,
                        SC::Comm,
                        SC,
                        SC::Domain,
                        D,
                    >(config, *other_domain, zeta, &mut circuit);

                    let first_point = circuit
                        .add_extension_constant(SC::Challenge::from(pcs.first_point(domain)));
                    let other_v_n =
                        vanishing_poly_at_point_circuit::<
                            InputProof,
                            OpeningProof,
                            SC::Comm,
                            SC,
                            SC::Domain,
                            D,
                        >(config, *other_domain, first_point, &mut circuit);
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
    InputProof: RecursiveVersion,
    OpeningProof: RecursiveVersion,
    Comm: RecursiveVersion,
    SC: RecursiveStarkGenerationConfig<InputProof, OpeningProof, D>,
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
    <SC as RecursiveStarkGenerationConfig<InputProof, OpeningProof, D>>::Pcs:
        PcsRecursiveVerif<InputProof, OpeningProof, Comm, Domain, SC::Val, SC::Challenge, D>,
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
