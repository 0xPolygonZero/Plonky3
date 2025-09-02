use p3_commit::{BatchOpening, Mmcs, PolynomialSpace, testing::TrivialPcs};
use p3_dft::TwoAdicSubgroupDft;
use p3_field::{
    ExtensionField, Field, TwoAdicField, coset::TwoAdicMultiplicativeCoset,
    extension::BinomiallyExtendable,
};
use p3_fri::TwoAdicFriPcs;
use p3_util::log2_strict_usize;

use crate::{
    circuit_builder::{
        CircuitBuilder, ExtensionWireId, WireId,
        gates::arith_gates::{MulExtensionGate, SubExtensionGate},
    },
    verifier::{
        circuit_verifier::ProofWires,
        recursive_traits::{
            CommitForRecursiveVerif, CommitRecursiveVerif, PcsRecursiveVerif, RecursiveLagrangeSels,
        },
    },
};

impl<Dft, InputMmcs: Mmcs<F>, FriMmcs, F: TwoAdicField, EF, const D: usize>
    PcsRecursiveVerif<
        Vec<BatchOpening<F, InputMmcs>>,
        InputMmcs::Commitment,
        TwoAdicMultiplicativeCoset<F>,
        F,
        EF,
        D,
    > for TwoAdicFriPcs<F, Dft, InputMmcs, FriMmcs>
where
    F: BinomiallyExtendable<D>,
    InputMmcs::Commitment: CommitRecursiveVerif + CommitForRecursiveVerif<F>,
    EF: ExtensionField<F>,
{
    fn get_challenges_circuit(
        circuit: &mut CircuitBuilder<F, D>,
        proof_wires: &ProofWires<D, InputMmcs::Commitment, Vec<BatchOpening<F, InputMmcs>>>,
    ) -> Vec<ExtensionWireId<D>> {
        let num_challenges = 1
            + proof_wires.fri_proof.commit_phase_commits.len()
            + proof_wires.fri_proof.query_proofs.len();

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
            &InputMmcs::Commitment,
            Vec<(
                TwoAdicMultiplicativeCoset<F>,
                Vec<([usize; D], Vec<[usize; D]>)>,
            )>,
        )],
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
            Vec<BatchOpening<F, InputMmcs>>,
            <InputMmcs as Mmcs<F>>::Commitment,
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

impl<F: Field> CommitForRecursiveVerif<F> for TrivialCommitment<F> {
    // fn get_commit_challenges_circuit(circuit: &mut CircuitBuilder<F, D>) -> Vec<WireId> {
    //     vec![]
    // }

    fn get_values(&self) -> Vec<F> {
        self.iter().flatten().cloned().collect()
    }
}

impl CommitRecursiveVerif for RecursiveTrivialCommitment {
    // fn get_commit_challenges_circuit(circuit: &mut CircuitBuilder<F, D>) -> Vec<WireId> {
    //     vec![]
    // }

    fn get_wires(&self) -> Vec<WireId> {
        self.iter().flatten().cloned().collect()
    }
}

impl<Dft, F: TwoAdicField, EF, const D: usize>
    PcsRecursiveVerif<(), RecursiveTrivialCommitment, TwoAdicMultiplicativeCoset<F>, F, EF, D>
    for TrivialPcs<F, Dft>
where
    Dft: TwoAdicSubgroupDft<F>,
    F: BinomiallyExtendable<D>,
    EF: ExtensionField<F>,
{
    fn get_challenges_circuit(
        circuit: &mut CircuitBuilder<F, D>,
        proof_wires: &ProofWires<D, RecursiveTrivialCommitment, ()>,
    ) -> Vec<ExtensionWireId<D>> {
        let num_challenges = 1
            + proof_wires.fri_proof.commit_phase_commits.len()
            + proof_wires.fri_proof.query_proofs.len();

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
            &RecursiveTrivialCommitment,
            Vec<(
                TwoAdicMultiplicativeCoset<F>,
                Vec<(ExtensionWireId<D>, Vec<ExtensionWireId<D>>)>,
            )>,
        )],
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
        let exp = circuit.add_extension_constant(EF::from_usize(
            <TrivialPcs<F, Dft> as PcsRecursiveVerif<
                (),
                RecursiveTrivialCommitment,
                TwoAdicMultiplicativeCoset<F>,
                F,
                EF,
                D,
            >>::size(self, domain),
        ));

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
