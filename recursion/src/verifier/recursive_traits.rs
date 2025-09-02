use p3_field::extension::BinomiallyExtendable;
use p3_field::{ExtensionField, Field};
use p3_symmetric::Hash;

use crate::circuit_builder::{ChallengeWireId, CircuitBuilder, WireId};
use crate::verifier::circuit_verifier::ProofWires;

pub trait CommitRecursiveVerif<F: Field, const D: usize> {
    fn get_commit_challenges_circuit(circuit_builder: &mut CircuitBuilder<F, D>) -> Vec<WireId>;
}

impl<F: Field, W, const DIGEST_ELEMS: usize, const D: usize> CommitRecursiveVerif<F, D>
    for Hash<F, W, DIGEST_ELEMS>
{
    fn get_commit_challenges_circuit(circuit_builder: &mut CircuitBuilder<F, D>) -> Vec<WireId> {
        (0..DIGEST_ELEMS)
            .map(|_| circuit_builder.new_wire())
            .collect()
    }
}

pub trait RecursiveStarkGenerationConfig<InputProof, const D: usize> {
    type Val: Field + BinomiallyExtendable<D>;
    type Domain: Copy;
    type Challenge: ExtensionField<Self::Val>;
    type Comm: CommitRecursiveVerif<Self::Val, D>;
    type Pcs: PcsRecursiveVerif<InputProof, Self::Comm, Self::Domain, Self::Val, Self::Challenge, D>;

    fn pcs(&self) -> Self::Pcs;

    fn is_zk(&self) -> usize;
}

pub trait PcsRecursiveVerif<
    InputProof,
    Comm: CommitRecursiveVerif<F, D>,
    Domain,
    F: Field,
    EF,
    const D: usize,
> where
    EF: ExtensionField<F>,
{
    fn get_challenges_circuit(
        circuit: &mut CircuitBuilder<F, D>,
        proof_wires: &ProofWires<F, D, Comm, InputProof>,
    ) -> Vec<ChallengeWireId<D>>;

    fn verify_circuit(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        challenges: &[ChallengeWireId<D>],
        commitments_with_opening_points: &[(
            &Comm,
            Vec<(Domain, Vec<([usize; D], Vec<[usize; D]>)>)>,
        )],
    );

    fn selectors_at_point_circuit(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        domain: &Domain,
        point: &ChallengeWireId<D>,
    ) -> RecursiveLagrangeSels<D>;

    fn natural_domain_for_degree(&self, degree: usize) -> Domain;

    fn create_disjoint_domain(&self, trace_domain: Domain, degree: usize) -> Domain;

    fn split_domains(&self, trace_domain: &Domain, degree: usize) -> Vec<Domain>;

    fn size(&self, trace_domain: &Domain) -> usize;

    fn first_point(&self, trace_domain: &Domain) -> F;
}

pub struct RecursiveLagrangeSels<const D: usize> {
    pub is_first_row: ChallengeWireId<D>,
    pub is_last_row: ChallengeWireId<D>,
    pub is_transition: ChallengeWireId<D>,
    pub inv_vanishing: ChallengeWireId<D>,
}

pub trait RecursiveAir<F: Field, const D: usize> {
    type Var: Clone;

    fn width(&self) -> usize;

    fn eval_folded_circuit(
        &self,
        builder: &mut CircuitBuilder<F, D>,
        sels: &RecursiveLagrangeSels<D>,
        alpha: &ChallengeWireId<D>,
        public_values: &[WireId],
    ) -> ChallengeWireId<D>;

    fn get_log_quotient_degree(
        &self,
        preprocessed_width: usize,
        num_public_values: usize,
        is_zk: usize,
    ) -> usize;
}
