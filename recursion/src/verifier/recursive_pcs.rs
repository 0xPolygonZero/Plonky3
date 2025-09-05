use std::{array, iter, marker::PhantomData};

use crate::{
    circuit_builder::{
        CircuitBuilder, ExtensionWireId, WireId,
        gates::arith_gates::{MulExtensionGate, SubExtensionGate},
    },
    verifier::recursive_traits::{PcsGeneration, Recursive, RecursiveLagrangeSels, RecursiveMmcs},
};
use p3_challenger::{CanObserve, CanSampleBits, FieldChallenger, GrindingChallenger};
use p3_commit::{BatchOpening, ExtensionMmcs, Mmcs, Pcs, PolynomialSpace};
use p3_field::{
    BasedVectorSpace, ExtensionField, Field, PackedValue, PrimeCharacteristicRing, TwoAdicField,
    coset::TwoAdicMultiplicativeCoset, extension::BinomiallyExtendable,
};
use p3_fri::{CommitPhaseProofStep, FriProof, QueryProof, TwoAdicFriPcs};
use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::Hash;
use p3_uni_stark::{StarkGenericConfig, Val};
use p3_util::log2_strict_usize;

// Define all the recursion types necessary for Fri and the TwoAdicFriPcs.
// Also, we implement `RecursionVersion` for these structures.

#[derive(Clone)]
pub struct FriProofWires<
    F: Field,
    EF: ExtensionField<F>,
    RecMmcs: RecursiveMmcs<F, D>,
    InputProof: Recursive<F, D>,
    Witness: Recursive<F, D>,
    const D: usize,
> {
    pub commit_phase_commits: Vec<RecMmcs::Commitment>,
    pub query_proofs: Vec<QueryProofWires<F, InputProof, RecMmcs, D>>,
    pub final_poly: Vec<ExtensionWireId<D>>,
    pub pow_witness: Witness,
    _phantom: PhantomData<EF>,
}

#[derive(Clone)]
pub struct QueryProofWires<
    F: Field,
    InputProof: Recursive<F, D>,
    RecMmcs: RecursiveMmcs<F, D>,
    const D: usize,
> {
    pub input_proof: InputProof,
    pub commit_phase_openings: Vec<CommitPhaseProofStepWires<F, RecMmcs, D>>,
}

#[derive(Clone)]
pub struct CommitPhaseProofStepWires<F: Field, RecMmcs: RecursiveMmcs<F, D>, const D: usize> {
    pub sibling_value: WireId,
    pub opening_proof: RecMmcs::Proof,
}

#[derive(Clone)]
pub struct BatchOpeningWires<F: Field, InnerProof: Recursive<F, D>, const D: usize> {
    /// The opened row values from each matrix in the batch.
    /// Each inner vector corresponds to one matrix.
    pub opened_values: Vec<Vec<WireId>>,
    /// The proof showing the values are valid openings.
    pub opening_proof: Vec<InnerProof>,
    _phantom: PhantomData<F>,
}

impl<F: Field, RecMmcs: RecursiveMmcs<F, D>, const D: usize> Recursive<F, D>
    for CommitPhaseProofStepWires<F, RecMmcs, D>
{
    type Input = CommitPhaseProofStep<F, RecMmcs::Input>;

    fn new(circuit: &mut CircuitBuilder<F, D>, lens: &mut impl Iterator<Item = usize>) -> Self {
        let sibling_value = circuit.new_wire();
        let opening_proof = <RecMmcs::Proof as Recursive<F, D>>::new(circuit, lens);
        Self {
            sibling_value,
            opening_proof,
        }
    }

    fn get_values(input: Self::Input) -> Vec<F> {
        let CommitPhaseProofStep {
            sibling_value,
            opening_proof,
        } = input;

        let mut values = vec![sibling_value];
        values.extend(<RecMmcs::Proof as Recursive<F, D>>::get_values(
            opening_proof,
        ));
        values
    }

    fn set_wires(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        input: Self::Input,
    ) -> Result<(), crate::circuit_builder::CircuitError> {
        let CommitPhaseProofStep {
            sibling_value,
            opening_proof,
        } = input;

        let CommitPhaseProofStepWires {
            sibling_value: sibling_value_wire,
            opening_proof: opening_proof_wires,
        } = self;

        circuit.set_wire_value(*sibling_value_wire, sibling_value)?;
        <RecMmcs::Proof as Recursive<F, D>>::set_wires(
            opening_proof_wires,
            circuit,
            opening_proof,
        )?;

        Ok(())
    }

    fn num_challenges(&self) -> usize {
        0
    }
}

impl<F: Field, InputProof: Recursive<F, D>, RecMmcs: RecursiveMmcs<F, D>, const D: usize>
    Recursive<F, D> for QueryProofWires<F, InputProof, RecMmcs, D>
{
    type Input = QueryProof<F, RecMmcs::Input, InputProof::Input>;

    fn new(circuit: &mut CircuitBuilder<F, D>, lens: &mut impl Iterator<Item = usize>) -> Self {
        let input_proof = InputProof::new(circuit, lens);
        let num_commit_phase_openings = lens.next().unwrap();
        let mut commit_phase_openings = Vec::with_capacity(num_commit_phase_openings);
        for _ in 0..num_commit_phase_openings {
            commit_phase_openings.push(CommitPhaseProofStepWires::<F, RecMmcs, D>::new(
                circuit, lens,
            ));
        }
        Self {
            input_proof,
            commit_phase_openings,
        }
    }

    fn get_values(input: Self::Input) -> Vec<F> {
        let QueryProof {
            input_proof,
            commit_phase_openings,
        } = input;

        let mut all_values = vec![];
        all_values.extend(<InputProof as Recursive<F, D>>::get_values(input_proof));
        all_values.extend(commit_phase_openings.iter().flat_map(|o| {
            <CommitPhaseProofStepWires<F, RecMmcs, D> as Recursive<F, D>>::get_values(o.clone())
        }));
        all_values
    }

    fn set_wires(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        input: Self::Input,
    ) -> Result<(), crate::circuit_builder::CircuitError> {
        let QueryProof {
            input_proof,
            commit_phase_openings,
        } = input;

        let QueryProofWires {
            input_proof: input_proof_wires,
            commit_phase_openings: commit_phase_openings_wires,
        } = self;

        <InputProof as Recursive<F, D>>::set_wires(&input_proof_wires, circuit, input_proof)?;
        for (cpo, w) in commit_phase_openings
            .iter()
            .zip(commit_phase_openings_wires.iter())
        {
            <CommitPhaseProofStepWires<F, RecMmcs, D> as Recursive<F, D>>::set_wires(
                w,
                circuit,
                cpo.clone(),
            )?;
        }

        Ok(())
    }

    fn num_challenges(&self) -> usize {
        0
    }
    // fn new(
    //     input_proof: InputProof,
    //     commit_phase_openings: Vec<CommitPhaseProofStepWires<RecursiveMmcsProof>>,
    // ) -> Self {
    //     Self {
    //         input_proof,
    //         commit_phase_openings,
    //         _phantom: PhantomData,
    //     }
    // }
}

pub type InputProofWires<F, Inner, const D: usize> = Vec<BatchOpeningWires<F, Inner, D>>;

type TwoadicFriProofWires<F, EF, RecMmcs, Inner, const D: usize> =
    FriProofWires<F, EF, RecMmcs, InputProofWires<F, Inner, D>, WireId, D>;

impl<
    F: Field,
    EF: ExtensionField<F>,
    RecMmcs: RecursiveMmcs<F, D>,
    InputProof: Recursive<F, D>,
    Witness: Recursive<F, D>,
    const D: usize,
> Recursive<F, D> for FriProofWires<F, EF, RecMmcs, InputProof, Witness, D>
{
    type Input = FriProof<F, RecMmcs::Input, Witness::Input, InputProof::Input>;

    fn get_values(input: Self::Input) -> Vec<F> {
        let FriProof {
            commit_phase_commits,
            query_proofs,
            final_poly,
            pow_witness,
        } = input;

        commit_phase_commits
            .iter()
            .flat_map(|c| {
                <<RecMmcs as RecursiveMmcs<F, D>>::Commitment as Recursive<F, D>>::get_values(
                    c.clone(),
                )
            })
            .chain(query_proofs.iter().flat_map(|c| {
                <QueryProofWires<F, InputProof, RecMmcs, D> as Recursive<F, D>>::get_values(
                    c.clone(),
                )
            }))
            .chain(final_poly.iter().cloned())
            .chain(<Witness as Recursive<F, D>>::get_values(pow_witness))
            .collect()
    }

    fn new(circuit: &mut CircuitBuilder<F, D>, lens: &mut impl Iterator<Item = usize>) -> Self {
        let num_commit_phase_commits = lens.next().unwrap();
        let mut commit_phase_commits = Vec::with_capacity(num_commit_phase_commits);
        for _ in 0..num_commit_phase_commits {
            commit_phase_commits.push(RecMmcs::Commitment::new(circuit, lens));
        }

        let num_query_proofs = lens.next().unwrap();
        let mut query_proofs = Vec::with_capacity(num_query_proofs);
        for _ in 0..num_query_proofs {
            query_proofs.push(QueryProofWires::<F, InputProof, RecMmcs, D>::new(
                circuit, lens,
            ));
        }
        // `lens` has been updated by the other structures. So the first element is indeed the length of the final polynomial.
        let final_poly_len = lens.next().unwrap();
        let mut final_poly = Vec::with_capacity(final_poly_len);
        for _ in 0..final_poly_len {
            final_poly.push(circuit.new_extension_wires());
        }
        Self {
            commit_phase_commits,
            query_proofs,
            final_poly,
            pow_witness: Witness::new(circuit, lens),
            _phantom: PhantomData,
        }
    }

    fn set_wires(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        input: Self::Input,
    ) -> Result<(), crate::circuit_builder::CircuitError> {
        let FriProof {
            commit_phase_commits,
            query_proofs,
            final_poly,
            pow_witness,
        } = input;

        let FriProofWires {
            commit_phase_commits: commit_phase_commits_wires,
            query_proofs: query_proofs_wires,
            final_poly: final_poly_wires,
            pow_witness: pow_witness_wires,
            ..
        } = self;

        for (c, w) in commit_phase_commits
            .iter()
            .zip(commit_phase_commits_wires.iter())
        {
            w.set_wires(circuit, c.clone())?;
        }
        for (q, w) in query_proofs.iter().zip(query_proofs_wires.iter()) {
            w.set_wires(circuit, q.clone())?;
        }
        for (f, w) in final_poly.iter().zip(final_poly_wires.iter()) {
            let f_ext: [F; D] = f.as_basis_coefficients_slice().try_into().unwrap();
            circuit.set_extension_wires(*w, &f_ext)?;
        }
        pow_witness_wires.set_wires(circuit, pow_witness)?;

        Ok(())
    }

    fn num_challenges(&self) -> usize {
        0
    }

    // fn lens(input: &Self::Input) -> impl Iterator<Item = usize> {
    //     let FriProof {
    //         commit_phase_commits,
    //         query_proofs,
    //         final_poly,
    //         pow_witness,
    //     } = input;

    //     let mut lens = vec![commit_phase_commits.len()];
    //     lens.extend(
    //         commit_phase_commits
    //             .iter()
    //             .flat_map(|c| RecMmcs::Comm::lens(c)),
    //     );
    //     lens.push(query_proofs.len());
    //     lens.extend(query_proofs.iter().flat_map(|q| q.lens()));
    //     lens.push(final_poly.len());
    //     lens.extend(*pow_witness.lens());

    //     lens.into_iter()
    // }
}

pub type InputProof<F, InputMmcs> = Vec<BatchOpening<F, InputMmcs>>;

// impl<F: Field, InputMmcs: Mmcs<F>> ForRecursiveVersion<F> for BatchOpening<F, InputMmcs>
// where
//     InputMmcs::Proof: ForRecursiveVersion<F>,
// {
//     type RecVerif = InputProofWires<InnerProof<F, InputMmcs>>;

//     fn get_values(&self) -> Vec<F> {
//         let mut all_values = vec![];
//         for row in &self.opened_values {
//             all_values.extend(row.iter().cloned());
//         }
//         all_values.extend(self.opening_proof.get_values());
//         all_values
//     }

//     fn lens(&self) -> Vec<usize> {
//         let mut lens = vec![];
//         lens.push(self.opened_values.len());
//         for row in &self.opened_values {
//             lens.push(row.len());
//         }
//         lens.extend(self.opening_proof.lens());
//         lens
//     }
// }

// impl<Inner: RecursiveVersion> RecursiveVersion for InputProofWires<Inner> {
//     fn get_wires(&self) -> Vec<WireId> {
//         self.iter()
//             .flat_map(|batch_opening| {
//                 batch_opening.opened_values.iter().flatten().cloned().chain(
//                     batch_opening
//                         .opening_proof
//                         .iter()
//                         .flat_map(|o| o.get_wires()),
//                 )
//             })
//             .collect()
//     }

//     fn new<F: Field, const K: usize>(
//         circuit: &mut CircuitBuilder<F, K>,
//         lens: &mut impl Iterator<Item = usize>,
//     ) -> Self {
//         let batch_openings_len = lens.next().unwrap();
//         let mut batch_openings = Vec::with_capacity(batch_openings_len);
//         for _ in 0..batch_openings_len {
//             let num_opened_rows = lens.next().unwrap();
//             let num_opened_cols = lens.next().unwrap();
//             let mut opened_values_rows = Vec::with_capacity(num_opened_rows);
//             for _ in 0..num_opened_rows {
//                 let mut opened_values_cols = Vec::with_capacity(num_opened_cols);
//                 for _ in 0..num_opened_cols {
//                     opened_values_cols.push(circuit.new_wire());
//                 }
//                 opened_values_rows.push(opened_values_cols);
//             }

//             let num_opening_proofs = lens.next().unwrap();
//             let mut opening_proof = Vec::with_capacity(num_opening_proofs);
//             for _ in 0..num_opening_proofs {
//                 opening_proof.push(Inner::new(circuit, lens));
//             }

//             batch_openings.push(BatchOpeningWires {
//                 opened_values: opened_values_rows,
//                 opening_proof,
//             })
//         }

//         batch_openings
//     }
// }

// type InnerProof<F, InputProof> = <<InputProof as ForRecursiveVersion<F>>::Proof as ForRecursiveVersion<F>>::RecVerif;

///////////// HERE ///////////////
// type InnerCommitment<F, InputMmcs> =
//     <<InputMmcs as Mmcs<F>>::Commitment as ForRecursiveVersion<F>>::RecVerif;

// type OpeningProof<F, InputMmcs> = FriProofWires<
//     InnerCommitment<F, InputMmcs>,
//     InnerProof<F, InputMmcs>,
//     InputProofWires<InnerProof<F, InputMmcs>>,
//     WireId,
// >;

// impl RecursiveVersion for WireId {
//     fn get_wires(&self) -> Vec<WireId> {
//         vec![*self]
//     }

//     fn new<F: Field, const D: usize>(
//         circuit: &mut CircuitBuilder<F, D>,
//         _lens: &mut impl Iterator<Item = usize>,
//     ) -> Self {
//         circuit.new_wire()
//     }
// }

// impl<SC: StarkGenericConfig, Dft, InputMmcs, FriMmcs>
//     PcsGeneration<
//         SC,
//         FriProof<SC::Challenge, FriMmcs, Val<SC>, Vec<BatchOpening<Val<SC>, InputMmcs>>>,
//     > for TwoAdicFriPcs<Val<SC>, Dft, InputMmcs, FriMmcs>
// where
//     SC::Challenge: ExtensionField<Val<SC>> + PrimeCharacteristicRing,
//     SC::Challenger: FieldChallenger<Val<SC>> + GrindingChallenger + CanObserve<FriMmcs::Commitment>,
//     FriMmcs: Mmcs<SC::Challenge>,
//     InputMmcs: Mmcs<Val<SC>>,
// {
//     fn generate_challenges<InputProof: ForRecursiveVersion<Val<SC>>, const D: usize>(
//         config: &SC,
//         challenger: &mut SC::Challenger,
//         coms_to_verify: &[(
//             <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment,
//             Vec<Vec<(SC::Challenge, Vec<SC::Challenge>)>>,
//         )],
//         opening_proof: &FriProof<
//             SC::Challenge,
//             FriMmcs,
//             Val<SC>,
//             Vec<BatchOpening<Val<SC>, InputMmcs>>,
//         >,
//     ) -> Vec<SC::Challenge> {
//         // Observe openings.
//         let mut challenges = vec![];
//         for (_, round) in coms_to_verify {
//             for mat in round {
//                 for (_, point) in mat {
//                     point.iter().for_each(|&opening| {
//                         challenger.observe_algebra_element(opening);
//                     })
//                 }
//             }
//         }

//         // Batch combination challenge `alpha`.
//         challenges.push(challenger.sample_algebra_element());

//         // Betas
//         let betas: Vec<SC::Challenge> = opening_proof
//             .commit_phase_commits
//             .iter()
//             .map(|comm| {
//                 // To match with the prover (and for security purposes),
//                 // we observe the commitment before sampling the challenge.
//                 challenger.observe(comm.clone());
//                 challenger.sample_algebra_element()
//             })
//             .collect();

//         challenges.extend(betas.iter());

//         // Observe all final polynomial cofficients.
//         opening_proof
//             .final_poly
//             .iter()
//             .for_each(|x| challenger.observe_algebra_element(*x));

//         let (log_blowup, log_final_poly_len) = config.pcs().get_log_blowup_final_height();
//         let log_max_height =
//             opening_proof.commit_phase_commits.len() + log_blowup + log_final_poly_len;

//         for _query_proof in &opening_proof.query_proofs {
//             challenges.push(SC::Challenge::from_usize(
//                 // Sample index. Note that `extra_query_bits` = 0 for the folding in this Pcs.
//                 challenger.sample_bits(log_max_height),
//             ));
//         }

//         challenges
//     }
// }

// // type ValMmcs<SC, Hash, Compress, const DIGEST_ELEMS: usize> = MerkleTreeMmcs<
// //     <Val<SC> as Field>::Packing,
// //     <Val<SC> as Field>::Packing,
// //     Hash,
// //     Compress,
// //     DIGEST_ELEMS,
// // >;
// // type ChallengeMmcs<SC, Hash, Compress, const DIGEST_ELEMS: usize> = ExtensionMmcs<
// //     Val<SC>,
// //     <SC as StarkGenericConfig>::Challenge,
// //     ValMmcs<SC, Hash, Compress, DIGEST_ELEMS>,
// // >;

type HashProofWires<const DIGEST_ELEMS: usize> = Vec<[WireId; DIGEST_ELEMS]>;
type ValMmcsProof<PW: PackedValue, const DIGEST_ELEMS: usize> = Vec<[PW::Value; DIGEST_ELEMS]>;

// impl<const DIGEST_ELEMS: usize> RecursiveVersion for HashProofWires<DIGEST_ELEMS> {
//     fn get_wires(&self) -> Vec<WireId> {
//         self.iter().flat_map(|h| h.to_vec()).collect()
//     }

//     fn new<F: Field, const D: usize>(
//         circuit: &mut CircuitBuilder<F, D>,
//         lens: &mut impl Iterator<Item = usize>,
//     ) -> Self {
//         let proof_len = lens.next().unwrap();
//         let mut proof = Vec::with_capacity(proof_len);
//         for _ in 0..proof_len {
//             proof.push(array::from_fn(|_| circuit.new_wire()));
//         }
//         proof
//     }
// }

// impl<F: Field, const DIGEST_ELEMS: usize> ForRecursiveVersion<F>
//     for ValMmcsProof<<F as Field>::Packing, DIGEST_ELEMS>
// {
//     type RecVerif = [WireId; DIGEST_ELEMS];

//     fn get_values(&self) -> Vec<F> {
//         self.iter().flat_map(|x| x.to_vec()).collect()
//     }

//     fn lens(&self) -> Vec<usize> {
//         vec![self.len()]
//     }
// }

type HashWires<const DIGEST_ELEMS: usize> = [WireId; DIGEST_ELEMS];
type ValMmcsCommitment<P: PackedValue, PW: PackedValue, const DIGEST_ELEMS: usize> =
    Hash<P::Value, PW::Value, DIGEST_ELEMS>;

// impl<const DIGEST_ELEMS: usize> RecursiveVersion for HashWires<DIGEST_ELEMS> {
//     fn get_wires(&self) -> Vec<WireId> {
//         self.to_vec()
//     }

//     fn new<F: Field, const D: usize>(
//         circuit: &mut CircuitBuilder<F, D>,
//         _lens: &mut impl Iterator<Item = usize>,
//     ) -> Self {
//         array::from_fn(|_| circuit.new_wire())
//     }
// }

// impl<F: Field, const DIGEST_ELEMS: usize> ForRecursiveVersion<F>
//     for ValMmcsCommitment<<F as Field>::Packing, <F as Field>::Packing, DIGEST_ELEMS>
// {
//     type RecVerif = [WireId; DIGEST_ELEMS];

//     fn get_values(&self) -> Vec<F> {
//         self.into_iter().collect::<Vec<_>>()
//     }

//     fn lens(&self) -> Vec<usize> {
//         vec![]
//     }
// }

// pub fn get_lens<SC: StarkGenericConfig, InputMmcs, FriMmcs>(
//     proof: &FriProof<SC::Challenge, FriMmcs, Val<SC>, Vec<BatchOpening<Val<SC>, InputMmcs>>>,
// ) -> Vec<usize>
// where
//     SC::Challenge: ExtensionField<Val<SC>> + PrimeCharacteristicRing,
//     SC::Challenger: FieldChallenger<Val<SC>> + GrindingChallenger + CanObserve<FriMmcs::Commitment>,
//     FriMmcs: Mmcs<SC::Challenge>,
//     InputMmcs: Mmcs<Val<SC>>,
//     InputMmcs::Commitment: ForRecursiveVersion<Val<SC>>,
//     FriMmcs::Commitment: ForRecursiveVersion<Val<SC>>,
//     InputMmcs::Proof: ForRecursiveVersion<Val<SC>>,
//     FriMmcs::Proof: ForRecursiveVersion<Val<SC>>,
// {
//     let FriProof {
//         commit_phase_commits,
//         query_proofs,
//         final_poly,
//         pow_witness: _,
//     } = proof;

//     let mut lens = vec![];
//     lens.push(commit_phase_commits.len());
//     for cpc in commit_phase_commits {
//         lens.extend(cpc.lens());
//     }

//     lens.push(query_proofs.len());
//     for qp in query_proofs {
//         lens.push(qp.commit_phase_openings.len());
//         for cpo in &qp.commit_phase_openings {
//             lens.extend(cpo.opening_proof.lens());
//         }
//         for bo in &qp.input_proof {
//             lens.push(bo.opened_values.len());
//             for ov in &bo.opened_values {
//                 lens.push(ov.len());
//             }
//             lens.extend(bo.opening_proof.lens());
//         }
//     }
//     lens.push(final_poly.len());
//     lens
// }

// // Implementing `PcsRecursiveVerif` for `TwoAdicFriPcs`.
// impl<Dft, InputMmcs: Mmcs<F>, FriMmcs, F: TwoAdicField, EF, const D: usize>
//     PcsRecursiveVerif<
//         InputProofWires<InnerProof<F, InputMmcs>>,
//         OpeningProof<F, InputMmcs>,
//         InnerCommitment<F, InputMmcs>,
//         TwoAdicMultiplicativeCoset<F>,
//         F,
//         EF,
//         D,
//     > for TwoAdicFriPcs<F, Dft, InputMmcs, FriMmcs>
// where
//     F: BinomiallyExtendable<D>,
//     InputMmcs::Commitment: ForRecursiveVersion<F>,
//     InputMmcs::Proof: ForRecursiveVersion<F>,
//     EF: ExtensionField<F>,
// {
//     type RecursiveProof = FriProofWires<
//         InnerCommitment<F, InputMmcs>,
//         InnerProof<F, InputMmcs>,
//         Vec<BatchOpeningWires<InnerProof<F, InputMmcs>>>,
//         WireId,
//     >;

//     fn get_challenges_circuit(
//         circuit: &mut CircuitBuilder<F, D>,
//         proof_wires: &ProofWires<D, InnerCommitment<F, InputMmcs>, OpeningProof<F, InputMmcs>>,
//     ) -> Vec<ExtensionWireId<D>> {
//         let num_challenges = 1
//             + proof_wires.opening_proof.commit_phase_commits.len()
//             + proof_wires.opening_proof.query_proofs.len();

//         let mut challenges = Vec::with_capacity(num_challenges);
//         for _ in 0..num_challenges {
//             challenges.push(circuit.new_extension_wires());
//         }

//         challenges
//     }

//     fn verify_circuit(
//         &self,
//         _circuit: &mut CircuitBuilder<F, D>,
//         _challenges: &[ExtensionWireId<D>],
//         _commitments_with_opening_points: &[(
//             &InnerCommitment<F, InputMmcs>,
//             Vec<(
//                 TwoAdicMultiplicativeCoset<F>,
//                 Vec<([usize; D], Vec<[usize; D]>)>,
//             )>,
//         )],
//         _opening_proof: &OpeningProof<F, InputMmcs>,
//     ) {
//         // For now, the verification doesn't do anything: we need the implementation of the FRI table first.
//         // Regarding the challenges, we only need interactions with the sponge tables.
//     }

//     fn selectors_at_point_circuit(
//         &self,
//         circuit: &mut CircuitBuilder<F, D>,
//         domain: &TwoAdicMultiplicativeCoset<F>,
//         point: &ExtensionWireId<D>,
//     ) -> RecursiveLagrangeSels<D> {
//         // Constants that we will need.
//         let shift_inv = circuit.add_extension_constant(EF::from(domain.shift_inverse()));
//         let one = circuit.add_extension_constant(EF::from(F::ONE));
//         let subgroup_gen_inv =
//             circuit.add_extension_constant(EF::from(domain.subgroup_generator().inverse()));
//         let exp = circuit.add_extension_constant(EF::from_usize(<TwoAdicFriPcs<
//             F,
//             Dft,
//             InputMmcs,
//             FriMmcs,
//         > as PcsRecursiveVerif<
//             InputProofWires<InnerProof<F, InputMmcs>>,
//             OpeningProof<F, InputMmcs>,
//             InnerCommitment<F, InputMmcs>,
//             TwoAdicMultiplicativeCoset<F>,
//             F,
//             EF,
//             D,
//         >>::size(self, domain)));

//         // Unshifted and z_h
//         let unshifted_point: [usize; D] = circuit.new_extension_wires();
//         MulExtensionGate::add_to_circuit(circuit, shift_inv, *point, unshifted_point);
//         let us_exp = circuit.new_extension_wires();
//         MulExtensionGate::add_to_circuit(circuit, unshifted_point, exp, us_exp);
//         let z_h = circuit.new_extension_wires();
//         SubExtensionGate::add_to_circuit(circuit, us_exp, one, z_h);

//         // Denominators
//         let us_minus_one = circuit.new_extension_wires();
//         SubExtensionGate::add_to_circuit(circuit, unshifted_point, one, us_minus_one);
//         let us_minus_gen_inv = circuit.new_extension_wires();
//         SubExtensionGate::add_to_circuit(
//             circuit,
//             unshifted_point,
//             subgroup_gen_inv,
//             us_minus_gen_inv,
//         );

//         // Selectors
//         let is_first_row = circuit.new_extension_wires();
//         MulExtensionGate::add_to_circuit(circuit, us_minus_one, is_first_row, z_h);
//         let is_last_row = circuit.new_extension_wires();
//         MulExtensionGate::add_to_circuit(circuit, us_minus_gen_inv, is_last_row, z_h);
//         let is_transition = us_minus_gen_inv;
//         let inv_vanishing = circuit.new_extension_wires();
//         MulExtensionGate::add_to_circuit(circuit, z_h, inv_vanishing, one);

//         RecursiveLagrangeSels {
//             is_first_row,
//             is_last_row,
//             is_transition,
//             inv_vanishing,
//         }
//     }

//     fn natural_domain_for_degree(&self, degree: usize) -> TwoAdicMultiplicativeCoset<F> {
//         TwoAdicMultiplicativeCoset::new(F::ONE, log2_strict_usize(degree)).unwrap()
//     }

//     fn create_disjoint_domain(
//         &self,
//         trace_domain: TwoAdicMultiplicativeCoset<F>,
//         degree: usize,
//     ) -> TwoAdicMultiplicativeCoset<F> {
//         trace_domain.create_disjoint_domain(degree)
//     }

//     fn split_domains(
//         &self,
//         trace_domain: &TwoAdicMultiplicativeCoset<F>,
//         degree: usize,
//     ) -> Vec<TwoAdicMultiplicativeCoset<F>> {
//         trace_domain.split_domains(degree)
//     }

//     fn size(&self, trace_domain: &TwoAdicMultiplicativeCoset<F>) -> usize {
//         trace_domain.size()
//     }

//     fn first_point(&self, trace_domain: &TwoAdicMultiplicativeCoset<F>) -> F {
//         trace_domain.first_point()
//     }
// }
