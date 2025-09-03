use itertools::Itertools;
use p3_commit::Pcs;
use p3_field::{ExtensionField, Field, extension::BinomiallyExtendable};
use p3_uni_stark::{StarkGenericConfig, Val};

use crate::{
    circuit_builder::{CircuitBuilder, CircuitError, ExtensionWireId, WireId},
    verifier::circuit_verifier::ProofWires,
};

/// Trait for the recursive commitments. It is used to extract the wires.
pub trait RecursiveVersion {
    /// Returns a vec of field elements representing one commitment.
    fn get_wires(&self) -> Vec<WireId>;
}

/// Trait for the non-recursive commitments. It is used to extract the values in the same format as the wires.
pub trait ForRecursiveVersion<F: Field> {
    type RecVerif: RecursiveVersion;
    /// Returns a vec of field elements representing one commitment.
    fn get_values(&self) -> Vec<F>;

    fn set_wires<SC: StarkGenericConfig, const D: usize>(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        recursive_version: Self::RecVerif,
    ) -> Result<(), CircuitError> {
        let values = self.get_values();
        let wires = recursive_version.get_wires();
        for (v, w) in values.iter().zip_eq(wires.iter()) {
            circuit.set_wire_value(*w, *v)?;
        }
        Ok(())
    }
}

// Note: might not be useful after all.
pub trait RecursiveStarkGenerationConfig<
    InputProof: RecursiveVersion,
    OpeningProof: RecursiveVersion,
    const D: usize,
>
{
    type Val: Field + BinomiallyExtendable<D>;
    type Domain: Copy;
    type Challenge: ExtensionField<Self::Val>;
    type Comm: RecursiveVersion;
    type Pcs: PcsRecursiveVerif<
            InputProof,
            OpeningProof,
            Self::Comm,
            Self::Domain,
            Self::Val,
            Self::Challenge,
            D,
        >;

    fn pcs(&self) -> Self::Pcs;

    fn is_zk(&self) -> usize;
}

/// Trait which defines the methods necessary
/// for a Pcs to generate values for associated wires.
pub trait PcsGeneration<SC: StarkGenericConfig, OpeningProof> {
    fn generate_challenges<InputProof: ForRecursiveVersion<Val<SC>>, const D: usize>(
        config: &SC,
        challenger: &mut SC::Challenger,
        coms_to_verify: &[(
            <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment,
            Vec<Vec<(SC::Challenge, Vec<SC::Challenge>)>>,
        )],
        opening_proof: &OpeningProof,
    ) -> Vec<SC::Challenge>;
}

/// Traits including the methods necessary for the recursive version of Pcs.
pub trait PcsRecursiveVerif<
    InputProof: RecursiveVersion,
    OpeningProof: RecursiveVersion,
    Comm: RecursiveVersion,
    Domain,
    F: Field,
    EF,
    const D: usize,
> where
    EF: ExtensionField<F>,
{
    type RecursiveProof;

    /// Creates new wires for all the challenges necessary when computing the Pcs.
    fn get_challenges_circuit(
        circuit: &mut CircuitBuilder<F, D>,
        proof_wires: &ProofWires<D, Comm, OpeningProof>,
    ) -> Vec<ExtensionWireId<D>>;

    /// Adds the circuit which verifies the Pcs computation.
    fn verify_circuit(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        challenges: &[ExtensionWireId<D>],
        commitments_with_opening_points: &[(
            &Comm,
            Vec<(Domain, Vec<([usize; D], Vec<[usize; D]>)>)>,
        )],
        opening_proof: &OpeningProof,
    );

    /// Computes wire selectors at `point` in the circuit.
    fn selectors_at_point_circuit(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        domain: &Domain,
        point: &ExtensionWireId<D>,
    ) -> RecursiveLagrangeSels<D>;

    /// Computes a domain given the degree. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    fn natural_domain_for_degree(&self, degree: usize) -> Domain;

    /// Computes a disjoint domain given the degree and the current domain. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    fn create_disjoint_domain(&self, trace_domain: Domain, degree: usize) -> Domain;

    /// Split a domain given the degree and the current domain. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    fn split_domains(&self, trace_domain: &Domain, degree: usize) -> Vec<Domain>;

    /// Returns the size of the domain. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    fn size(&self, trace_domain: &Domain) -> usize;

    /// Returns the first point in the domain. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    fn first_point(&self, trace_domain: &Domain) -> F;
}

/// Circuit version of the `LangrangeSelectors`.
pub struct RecursiveLagrangeSels<const D: usize> {
    pub is_first_row: ExtensionWireId<D>,
    pub is_last_row: ExtensionWireId<D>,
    pub is_transition: ExtensionWireId<D>,
    pub inv_vanishing: ExtensionWireId<D>,
}

/// Trait including methods necessary to compute the verification of an AIR's constraints,
/// as well as AIR-specific methods used in the full verification circuit.
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
