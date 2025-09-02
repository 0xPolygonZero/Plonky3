use p3_field::{ExtensionField, Field, extension::BinomiallyExtendable};

use crate::{
    circuit_builder::{ExtensionWireId, CircuitBuilder, WireId},
    verifier::circuit_verifier::ProofWires,
};

pub trait CommitRecursiveVerif {
    /// Returns a vec of field elements representing one commitment.
    fn get_wires(&self) -> Vec<WireId>;
}

pub trait CommitForRecursiveVerif<F: Field> {
    // fn get_commit_challenges_circuit(circuit: &mut CircuitBuilder<F, D>) -> Vec<WireId>;

    /// Returns a vec of field elements representing one commitment.
    fn get_values(&self) -> Vec<F>;
}

pub trait RecursiveStarkGenerationConfig<InputProof, const D: usize> {
    type Val: Field + BinomiallyExtendable<D>;
    type Domain: Copy;
    type Challenge: ExtensionField<Self::Val>;
    type Comm: CommitRecursiveVerif;
    type Pcs: PcsRecursiveVerif<InputProof, Self::Comm, Self::Domain, Self::Val, Self::Challenge, D>;

    fn pcs(&self) -> Self::Pcs;

    fn is_zk(&self) -> usize;
}

pub trait PcsRecursiveVerif<
    InputProof,
    Comm: CommitRecursiveVerif,
    Domain,
    F: Field,
    EF,
    const D: usize,
> where
    EF: ExtensionField<F>,
{
    fn get_challenges_circuit(
        circuit: &mut CircuitBuilder<F, D>,
        proof_wires: &ProofWires<D, Comm, InputProof>,
    ) -> Vec<ExtensionWireId<D>>;

    fn verify_circuit(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        challenges: &[ExtensionWireId<D>],
        commitments_with_opening_points: &[(
            &Comm,
            Vec<(Domain, Vec<([usize; D], Vec<[usize; D]>)>)>,
        )],
    );

    fn selectors_at_point_circuit(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        domain: &Domain,
        point: &ExtensionWireId<D>,
    ) -> RecursiveLagrangeSels<D>;

    fn natural_domain_for_degree(&self, degree: usize) -> Domain;

    fn create_disjoint_domain(&self, trace_domain: Domain, degree: usize) -> Domain;

    fn split_domains(&self, trace_domain: &Domain, degree: usize) -> Vec<Domain>;

    fn size(&self, trace_domain: &Domain) -> usize;

    fn first_point(&self, trace_domain: &Domain) -> F;
}

pub struct RecursiveLagrangeSels<const D: usize> {
    pub is_first_row: ExtensionWireId<D>,
    pub is_last_row: ExtensionWireId<D>,
    pub is_transition: ExtensionWireId<D>,
    pub inv_vanishing: ExtensionWireId<D>,
}

pub trait RecursiveAir<F: Field, const D: usize> {
    fn width(&self) -> usize;

    fn eval_folded_circuit<EF: ExtensionField<F>>(
        &self,
        builder: &mut CircuitBuilder<F, D>,
        sels: &RecursiveLagrangeSels<D>,
        alpha: &ExtensionWireId<D>,
        local_prep_values: &[ExtensionWireId<D>],
        next_prep_values: &[ExtensionWireId<D>],
        local_values: &[ExtensionWireId<D>],
        next_values: &[ExtensionWireId<D>],
        public_values: &[WireId],
    ) -> ExtensionWireId<D>;

    fn get_log_quotient_degree(
        &self,
        preprocessed_width: usize,
        num_public_values: usize,
        is_zk: usize,
    ) -> usize;
}
