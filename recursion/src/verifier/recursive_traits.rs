use p3_field::{ExtensionField, Field};
use p3_uni_stark::{Domain, StarkGenericConfig};

use crate::{
    circuit_builder::{ChallengeWireId, CircuitBuilder, WireId},
    verifier::circuit_verifier::ProofWires,
};

pub trait CommitRecursiveVerif {
    fn get_commit_challenges_circuit() -> Vec<WireId>;
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
    ) -> Vec<ChallengeWireId<D>>;

    fn verify_circuit(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        zeta: ChallengeWireId<D>,
        zeta_next: ChallengeWireId<D>,
        challenges: &[ChallengeWireId<D>],
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

pub trait RecursiveAir<SC: StarkGenericConfig, F: Field, const D: usize> {
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

pub trait PcsRecursiveGeneration<
    Challenger,
    InputProof,
    Comm: CommitRecursiveVerif,
    Domain,
    F: Field,
    EF,
    const D: usize,
> where
    EF: ExtensionField<F>,
{
    fn generate_challenges_circuit(
        circuit: &mut CircuitBuilder<F, D>,
        challenger: &mut Challenger,
        proof_wires: &ProofWires<D, Comm, InputProof>,
    ) -> Vec<ChallengeWireId<D>>;

    fn generate(&self, circuit: CircuitBuilder<F, D>, wires: &[ChallengeWireId<D>], inputs: &[EF]);
    fn generate_proof(
        &self,
        circuit: CircuitBuilder<F, D>,
        wires: &[ChallengeWireId<D>],
        inputs: &[EF],
    ) {
    }
}
