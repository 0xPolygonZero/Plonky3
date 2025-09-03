use p3_challenger::{CanObserve, CanSampleBits, FieldChallenger, GrindingChallenger};
use p3_commit::{BatchOpening, Mmcs, Pcs, PolynomialSpace};
use p3_field::{
    ExtensionField, Field, PrimeCharacteristicRing, TwoAdicField,
    coset::TwoAdicMultiplicativeCoset, extension::BinomiallyExtendable,
};
use p3_fri::{FriProof, TwoAdicFriPcs};
use p3_uni_stark::{StarkGenericConfig, Val};
use p3_util::log2_strict_usize;

use crate::{
    circuit_builder::{
        CircuitBuilder, ExtensionWireId, WireId,
        gates::arith_gates::{MulExtensionGate, SubExtensionGate},
    },
    verifier::{
        circuit_verifier::ProofWires,
        recursive_traits::{
            ForRecursiveVersion, PcsGeneration, PcsRecursiveVerif, RecursiveLagrangeSels,
            RecursiveVersion,
        },
    },
};

// Define all the recursion types necessary for Fri and the TwoAdicFriPcs.
// Also, we implement `RecursionVersion` for these structures.

#[derive(Clone)]
pub struct FriProofWires<
    const D: usize,
    Comm: RecursiveVersion,
    RecursiveMmcsProof: RecursiveVersion,
    InputProof: RecursiveVersion,
    Witness: RecursiveVersion,
> {
    pub commit_phase_commits: Vec<Comm>,
    pub query_proofs: Vec<QueryProofWires<InputProof, RecursiveMmcsProof>>,
    pub final_poly: Vec<WireId>,
    pub pow_witness: Witness,
}

impl<
    const D: usize,
    Comm: RecursiveVersion,
    RecursiveMmcsProof: RecursiveVersion,
    InputProof: RecursiveVersion,
    Witness: RecursiveVersion,
> RecursiveVersion for FriProofWires<D, Comm, RecursiveMmcsProof, InputProof, Witness>
{
    fn get_wires(&self) -> Vec<WireId> {
        let commit_wires = self.commit_phase_commits.iter().flat_map(|c| c.get_wires());
        let query_proofs = self.query_proofs.iter().flat_map(|q| q.get_wires());
        let final_poly = self.final_poly.iter().cloned();
        let pow_witness = self.pow_witness.get_wires();

        commit_wires
            .chain(query_proofs)
            .chain(final_poly)
            .chain(pow_witness)
            .collect()
    }
}

#[derive(Clone)]
pub struct QueryProofWires<InputProof: RecursiveVersion, RecursiveMmcsProof: RecursiveVersion> {
    pub input_proof: InputProof,
    pub commit_phase_openings: Vec<CommitPhaseProofStepWires<RecursiveMmcsProof>>,
}

impl<InputProof: RecursiveVersion, RecursiveMmcsProof: RecursiveVersion> RecursiveVersion
    for QueryProofWires<InputProof, RecursiveMmcsProof>
{
    fn get_wires(&self) -> Vec<WireId> {
        let input_proof_wires = self.input_proof.get_wires();
        let commit_wires = self
            .commit_phase_openings
            .iter()
            .map(|o| o.get_wires())
            .flatten();
        input_proof_wires.into_iter().chain(commit_wires).collect()
    }
}

#[derive(Clone)]
pub struct CommitPhaseProofStepWires<RecursiveMmcsProof> {
    pub sibling_value: WireId,
    pub opening_proof: Vec<RecursiveMmcsProof>,
}

impl<RecursiveMmcsProof: RecursiveVersion> RecursiveVersion
    for CommitPhaseProofStepWires<RecursiveMmcsProof>
{
    fn get_wires(&self) -> Vec<WireId> {
        std::iter::once(self.sibling_value)
            .chain(self.opening_proof.iter().flat_map(|p| p.get_wires()))
            .collect()
    }
}

#[derive(Clone)]
pub struct BatchOpeningWires<RecursiveMmcsProof> {
    /// The opened row values from each matrix in the batch.
    /// Each inner vector corresponds to one matrix.
    pub opened_values: Vec<Vec<WireId>>,
    /// The proof showing the values are valid openings.
    pub opening_proof: Vec<RecursiveMmcsProof>,
}

pub type InputProof<Inner> = Vec<BatchOpeningWires<Inner>>;

impl<Inner: RecursiveVersion> RecursiveVersion for InputProof<Inner> {
    fn get_wires(&self) -> Vec<WireId> {
        self.iter()
            .flat_map(|batch_opening| {
                batch_opening.opened_values.iter().flatten().cloned().chain(
                    batch_opening
                        .opening_proof
                        .iter()
                        .flat_map(|o| o.get_wires()),
                )
            })
            .collect()
    }
}

type InnerProof<F, InputMmcs> = <<InputMmcs as Mmcs<F>>::Proof as ForRecursiveVersion<F>>::RecVerif;

type InnerCommitment<F, InputMmcs> =
    <<InputMmcs as Mmcs<F>>::Commitment as ForRecursiveVersion<F>>::RecVerif;

type OpeningProof<const D: usize, F, InputMmcs> = FriProofWires<
    D,
    InnerCommitment<F, InputMmcs>,
    InnerProof<F, InputMmcs>,
    InputProof<InnerProof<F, InputMmcs>>,
    WireId,
>;

impl RecursiveVersion for WireId {
    fn get_wires(&self) -> Vec<WireId> {
        vec![*self]
    }
}

impl<SC: StarkGenericConfig, Dft, InputMmcs, FriMmcs>
    PcsGeneration<
        SC,
        FriProof<SC::Challenge, FriMmcs, Val<SC>, Vec<BatchOpening<Val<SC>, InputMmcs>>>,
    > for TwoAdicFriPcs<Val<SC>, Dft, InputMmcs, FriMmcs>
where
    SC::Challenge: ExtensionField<Val<SC>> + PrimeCharacteristicRing,
    SC::Challenger: FieldChallenger<Val<SC>> + GrindingChallenger + CanObserve<FriMmcs::Commitment>,
    FriMmcs: Mmcs<SC::Challenge>,
    InputMmcs: Mmcs<Val<SC>>,
{
    fn generate_challenges<InputProof: ForRecursiveVersion<Val<SC>>, const D: usize>(
        config: &SC,
        challenger: &mut SC::Challenger,
        coms_to_verify: &[(
            <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment,
            Vec<Vec<(SC::Challenge, Vec<SC::Challenge>)>>,
        )],
        opening_proof: &FriProof<
            SC::Challenge,
            FriMmcs,
            Val<SC>,
            Vec<BatchOpening<Val<SC>, InputMmcs>>,
        >,
    ) -> Vec<SC::Challenge> {
        // Observe openings.
        let mut challenges = vec![];
        for (_, round) in coms_to_verify {
            for mat in round {
                for (_, point) in mat {
                    point.iter().for_each(|&opening| {
                        challenger.observe_algebra_element(opening);
                    })
                }
            }
        }

        // Batch combination challenge `alpha`.
        challenges.push(challenger.sample_algebra_element());

        // Betas
        let betas: Vec<SC::Challenge> = opening_proof
            .commit_phase_commits
            .iter()
            .map(|comm| {
                // To match with the prover (and for security purposes),
                // we observe the commitment before sampling the challenge.
                challenger.observe(comm.clone());
                challenger.sample_algebra_element()
            })
            .collect();

        challenges.extend(betas.iter());

        // Observe all final polynomial cofficients.
        opening_proof
            .final_poly
            .iter()
            .for_each(|x| challenger.observe_algebra_element(*x));

        let (log_blowup, log_final_poly_len) = config.pcs().get_log_blowup_final_height();
        let log_max_height =
            opening_proof.commit_phase_commits.len() + log_blowup + log_final_poly_len;

        for _query_proof in &opening_proof.query_proofs {
            challenges.push(SC::Challenge::from_usize(
                // Sample index. Note that `extra_query_bits` = 0 for the folding in this Pcs.
                challenger.sample_bits(log_max_height),
            ));
        }

        challenges
    }
}

// Implementing `PcsRecursiveVerif` for `TwoAdicFriPcs`.
impl<Dft, InputMmcs: Mmcs<F>, FriMmcs, F: TwoAdicField, EF, const D: usize>
    PcsRecursiveVerif<
        InputProof<InnerProof<F, InputMmcs>>,
        OpeningProof<D, F, InputMmcs>,
        InnerCommitment<F, InputMmcs>,
        TwoAdicMultiplicativeCoset<F>,
        F,
        EF,
        D,
    > for TwoAdicFriPcs<F, Dft, InputMmcs, FriMmcs>
where
    F: BinomiallyExtendable<D>,
    InputMmcs::Commitment: ForRecursiveVersion<F>,
    InputMmcs::Proof: ForRecursiveVersion<F>,
    EF: ExtensionField<F>,
{
    type RecursiveProof = FriProofWires<
        D,
        InnerCommitment<F, InputMmcs>,
        InnerProof<F, InputMmcs>,
        Vec<BatchOpeningWires<InnerProof<F, InputMmcs>>>,
        WireId,
    >;
    fn get_challenges_circuit(
        circuit: &mut CircuitBuilder<F, D>,
        proof_wires: &ProofWires<D, InnerCommitment<F, InputMmcs>, OpeningProof<D, F, InputMmcs>>,
    ) -> Vec<ExtensionWireId<D>> {
        let num_challenges = 1
            + proof_wires.opening_proof.commit_phase_commits.len()
            + proof_wires.opening_proof.query_proofs.len();

        let mut challenges = Vec::with_capacity(num_challenges);
        for _ in 0..num_challenges {
            challenges.push(circuit.new_extension_wires());
        }

        challenges
    }

    fn verify_circuit(
        &self,
        _circuit: &mut CircuitBuilder<F, D>,
        _challenges: &[ExtensionWireId<D>],
        _commitments_with_opening_points: &[(
            &InnerCommitment<F, InputMmcs>,
            Vec<(
                TwoAdicMultiplicativeCoset<F>,
                Vec<([usize; D], Vec<[usize; D]>)>,
            )>,
        )],
        _opening_proof: &OpeningProof<D, F, InputMmcs>,
    ) {
        // For now, the verification doesn't do anything: we need the implementation of the FRI table first.
        // Regarding the challenges, we only need interactions with the sponge tables.
    }

    fn selectors_at_point_circuit(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        domain: &TwoAdicMultiplicativeCoset<F>,
        point: &ExtensionWireId<D>,
    ) -> RecursiveLagrangeSels<D> {
        // Constants that we will need.
        let shift_inv = circuit.add_extension_constant(EF::from(domain.shift_inverse()));
        let one = circuit.add_extension_constant(EF::from(F::ONE));
        let subgroup_gen_inv =
            circuit.add_extension_constant(EF::from(domain.subgroup_generator().inverse()));
        let exp = circuit.add_extension_constant(EF::from_usize(<TwoAdicFriPcs<
            F,
            Dft,
            InputMmcs,
            FriMmcs,
        > as PcsRecursiveVerif<
            InputProof<InnerProof<F, InputMmcs>>,
            OpeningProof<D, F, InputMmcs>,
            InnerCommitment<F, InputMmcs>,
            TwoAdicMultiplicativeCoset<F>,
            F,
            EF,
            D,
        >>::size(self, domain)));

        // Unshifted and z_h
        let unshifted_point: [usize; D] = circuit.new_extension_wires();
        MulExtensionGate::add_to_circuit(circuit, shift_inv, *point, unshifted_point);
        let us_exp = circuit.new_extension_wires();
        MulExtensionGate::add_to_circuit(circuit, unshifted_point, exp, us_exp);
        let z_h = circuit.new_extension_wires();
        SubExtensionGate::add_to_circuit(circuit, us_exp, one, z_h);

        // Denominators
        let us_minus_one = circuit.new_extension_wires();
        SubExtensionGate::add_to_circuit(circuit, unshifted_point, one, us_minus_one);
        let us_minus_gen_inv = circuit.new_extension_wires();
        SubExtensionGate::add_to_circuit(
            circuit,
            unshifted_point,
            subgroup_gen_inv,
            us_minus_gen_inv,
        );

        // Selectors
        let is_first_row = circuit.new_extension_wires();
        MulExtensionGate::add_to_circuit(circuit, us_minus_one, is_first_row, z_h);
        let is_last_row = circuit.new_extension_wires();
        MulExtensionGate::add_to_circuit(circuit, us_minus_gen_inv, is_last_row, z_h);
        let is_transition = us_minus_gen_inv;
        let inv_vanishing = circuit.new_extension_wires();
        MulExtensionGate::add_to_circuit(circuit, z_h, inv_vanishing, one);

        RecursiveLagrangeSels {
            is_first_row,
            is_last_row,
            is_transition,
            inv_vanishing,
        }
    }

    fn natural_domain_for_degree(&self, degree: usize) -> TwoAdicMultiplicativeCoset<F> {
        TwoAdicMultiplicativeCoset::new(F::ONE, log2_strict_usize(degree)).unwrap()
    }

    fn create_disjoint_domain(
        &self,
        trace_domain: TwoAdicMultiplicativeCoset<F>,
        degree: usize,
    ) -> TwoAdicMultiplicativeCoset<F> {
        trace_domain.create_disjoint_domain(degree)
    }

    fn split_domains(
        &self,
        trace_domain: &TwoAdicMultiplicativeCoset<F>,
        degree: usize,
    ) -> Vec<TwoAdicMultiplicativeCoset<F>> {
        trace_domain.split_domains(degree)
    }

    fn size(&self, trace_domain: &TwoAdicMultiplicativeCoset<F>) -> usize {
        trace_domain.size()
    }

    fn first_point(&self, trace_domain: &TwoAdicMultiplicativeCoset<F>) -> F {
        trace_domain.first_point()
    }
}

type TrivialCommitment<F> = Vec<Vec<F>>;
type RecursiveTrivialCommitment = Vec<Vec<WireId>>;

impl<F: Field> ForRecursiveVersion<F> for TrivialCommitment<F> {
    type RecVerif = RecursiveTrivialCommitment;

    fn get_values(&self) -> Vec<F> {
        self.iter().flatten().cloned().collect()
    }
}

impl RecursiveVersion for RecursiveTrivialCommitment {
    fn get_wires(&self) -> Vec<WireId> {
        self.iter().flatten().cloned().collect()
    }
}
