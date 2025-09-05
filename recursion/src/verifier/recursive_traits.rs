use std::marker::PhantomData;

use itertools::Itertools;
use p3_commit::{Mmcs, Pcs};
use p3_field::{BasedVectorSpace, ExtensionField, Field};
use p3_uni_stark::{Commitments, OpenedValues, Proof, StarkGenericConfig, Val};

use crate::circuit_builder::{CircuitBuilder, CircuitError, ExtensionWireId, WireId};

/// Structure representing all the wires necessary for an input proof.
#[derive(Clone)]
pub struct ProofWires<
    SC: StarkGenericConfig,
    Comm: Recursive<Val<SC>, D>,
    OpeningProof: Recursive<Val<SC>, D>,
    const D: usize,
> {
    pub commitments_wires: CommitmentWires<Val<SC>, Comm, D>,
    pub opened_values_wires: OpenedValuesWires<SC, D>,
    pub opening_proof: OpeningProof,
    pub degree_bits: usize,
}

#[derive(Clone)]
pub struct CommitmentWires<F: Field, Comm: Recursive<F, D>, const D: usize> {
    pub trace_wires: Comm,
    pub quotient_chunks_wires: Comm,
    pub random_commit: Option<Comm>,
    pub _phantom: PhantomData<F>,
}

// TODO: Move these structures to their respective crates.
#[derive(Clone)]
pub struct OpenedValuesWires<SC: StarkGenericConfig, const D: usize> {
    pub trace_local_wires: Vec<ExtensionWireId<D>>,
    pub trace_next_wires: Vec<ExtensionWireId<D>>,
    pub quotient_chunks_wires: Vec<Vec<ExtensionWireId<D>>>,
    pub random_wires: Option<Vec<ExtensionWireId<D>>>,
    _phantom: PhantomData<SC>,
}

pub trait Recursive<F: Field, const D: usize> {
    /// The nonrecursive type associated with the recursive type implementing the trait.
    type Input: Clone;

    /// Creates a new instance of the recursive type. `lens` corresponds to all the vector lengths necessary to build the structure. TODO: They can actually be deduced from StarkGenericConfig and `degree_bits`.
    fn new(circuit: &mut CircuitBuilder<F, D>, lens: &mut impl Iterator<Item = usize>) -> Self;

    /// Returns a vec of field elements representing the elements of the Input. Used (at least for now) to generate challenges.
    fn get_values(input: Self::Input) -> Vec<F>;

    /// Given an `Input` instance, sets the wires in the current instance to the values in `input`.
    fn set_wires(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        input: Self::Input,
    ) -> Result<(), CircuitError>;

    /// Returns the number of challenges necessary
    fn num_challenges(&self) -> usize;

    /// Creates new wires for all the necessary challenges.
    fn get_challenges(&self, circuit: &mut CircuitBuilder<F, D>) -> Vec<ExtensionWireId<D>> {
        let num_challenges = self.num_challenges();

        let mut challenges = Vec::with_capacity(num_challenges);
        for _ in 0..num_challenges {
            challenges.push(circuit.new_extension_wires());
        }

        challenges
    }

    // Temporary method used for testing for now. This should be changed into something more generic which relies as little as possible on the actual proof.
    // fn lens(input: &Self::Input) -> impl Iterator<Item = usize>;
}
// Note: might not be useful after all.
pub trait RecursiveStarkGenericConfig<
    SC: StarkGenericConfig,
    InputProof: Recursive<Val<SC>, D>,
    OpeningProof: Recursive<Val<SC>, D>,
    const D: usize,
>
{
    type Domain: Copy;
    type Comm: Recursive<Val<SC>, D>;
    type Pcs: RecursivePcs<SC, InputProof, OpeningProof, Self::Comm, Self::Domain, D>;

    fn pcs(&self) -> Self::Pcs;

    fn is_zk(&self) -> usize;
}

/// Trait which defines the methods necessary
/// for a Pcs to generate values for associated wires.
/// Generalize
pub trait PcsGeneration<SC: StarkGenericConfig, OpeningProof> {
    fn generate_challenges<InputProof: Recursive<Val<SC>, D>, const D: usize>(
        config: &SC,
        challenger: &mut SC::Challenger,
        coms_to_verify: &[(
            <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment,
            Vec<Vec<(SC::Challenge, Vec<SC::Challenge>)>>,
        )],
        opening_proof: &OpeningProof,
    ) -> Vec<SC::Challenge>;
}

// Can we not have this extra trait? Another possibility would be to have the `Input` in `Recursive` for `Commitment` be the entire Mmcs...
pub trait RecursiveMmcs<F: Field, const D: usize> {
    type Input: Mmcs<F>;
    type Commitment: Recursive<F, D, Input = <Self::Input as Mmcs<F>>::Commitment> + Clone;
    type Proof: Recursive<F, D, Input = <Self::Input as Mmcs<F>>::Proof> + Clone;
}

/// Traits including the methods necessary for the recursive version of Pcs.
/// Prepend Recursive
pub trait RecursivePcs<
    SC: StarkGenericConfig,
    InputProof: Recursive<Val<SC>, D>,
    OpeningProof: Recursive<Val<SC>, D>,
    Comm: Recursive<Val<SC>, D>,
    Domain,
    const D: usize,
>
{
    type RecursiveProof;

    /// Creates new wires for all the challenges necessary when computing the Pcs.
    fn get_challenges_circuit(
        circuit: &mut CircuitBuilder<Val<SC>, D>,
        proof_wires: &ProofWires<SC, Comm, OpeningProof, D>,
    ) -> Vec<ExtensionWireId<D>>;

    /// Adds the circuit which verifies the Pcs computation.
    fn verify_circuit(
        &self,
        circuit: &mut CircuitBuilder<Val<SC>, D>,
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
        circuit: &mut CircuitBuilder<Val<SC>, D>,
        domain: &Domain,
        point: &ExtensionWireId<D>,
    ) -> RecursiveLagrangeSels<D>;

    // /// Computes a domain given the degree. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    // fn natural_domain_for_degree(&self, degree: usize) -> Domain;

    /// Computes a disjoint domain given the degree and the current domain. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    fn create_disjoint_domain(&self, trace_domain: Domain, degree: usize) -> Domain;

    /// Split a domain given the degree and the current domain. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    fn split_domains(&self, trace_domain: &Domain, degree: usize) -> Vec<Domain>;

    /// Returns the size of the domain. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    fn size(&self, trace_domain: &Domain) -> usize;

    /// Returns the first point in the domain. This is the same as the original method in Pcs, but is also used in the verifier circuit.
    fn first_point(&self, trace_domain: &Domain) -> Val<SC>;
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

impl<
    SC: StarkGenericConfig + Clone,
    const D: usize,
    Comm: Recursive<Val<SC>, D, Input = <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment>,
    OpeningProof: Recursive<Val<SC>, D, Input = <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Proof>,
> Recursive<Val<SC>, D> for ProofWires<SC, Comm, OpeningProof, D>
{
    type Input = Proof<SC>;

    fn new(
        circuit: &mut CircuitBuilder<Val<SC>, D>,
        lens: &mut impl Iterator<Item = usize>,
    ) -> Self {
        let commitments_wires = CommitmentWires::new(circuit, lens);
        let opened_values_wires = OpenedValuesWires::new(circuit, lens);
        let opening_proof = OpeningProof::new(circuit, lens);
        let degree_bits = lens.next().unwrap();

        Self {
            commitments_wires,
            opened_values_wires,
            opening_proof,
            degree_bits,
        }
    }

    fn get_values(input: Self::Input) -> Vec<Val<SC>> {
        let Proof {
            commitments,
            opened_values,
            opening_proof,
            degree_bits: _,
        } = input;
        let mut values = vec![];
        values.extend::<Vec<Val<SC>>>(CommitmentWires::<Val<SC>, Comm, D>::get_values(commitments));
        values.extend(OpenedValuesWires::<SC, D>::get_values(opened_values));
        values.extend(OpeningProof::get_values(opening_proof));
        values
    }

    fn set_wires(
        &self,
        circuit: &mut CircuitBuilder<Val<SC>, D>,
        input: Self::Input,
    ) -> Result<(), CircuitError> {
        let Proof {
            commitments,
            opened_values,
            opening_proof,
            degree_bits,
        } = input;
        if degree_bits != self.degree_bits {
            return Err(CircuitError::DegreeBitsMismatch);
        }
        self.commitments_wires.set_wires(circuit, commitments)?;
        self.opened_values_wires.set_wires(circuit, opened_values)?;
        self.opening_proof.set_wires(circuit, opening_proof)?;
        Ok(())
    }

    fn num_challenges(&self) -> usize {
        self.commitments_wires.num_challenges()
            + self.opened_values_wires.num_challenges()
            + self.opening_proof.num_challenges()
    }
}

impl<F: Field, const D: usize, Comm> Recursive<F, D> for CommitmentWires<F, Comm, D>
where
    Comm: Recursive<F, D>,
{
    type Input = Commitments<Comm::Input>;

    fn new(circuit: &mut CircuitBuilder<F, D>, lens: &mut impl Iterator<Item = usize>) -> Self {
        let trace_wires = Comm::new(circuit, lens);
        let quotient_chunks_wires = Comm::new(circuit, lens);
        let random_commit = if lens.next().is_some() {
            Some(Comm::new(circuit, lens))
        } else {
            None
        };
        Self {
            trace_wires,
            quotient_chunks_wires,
            random_commit,
            _phantom: PhantomData,
        }
    }

    fn get_values(input: Self::Input) -> Vec<F> {
        let Commitments {
            trace,
            quotient_chunks,
            random,
        } = input;

        let mut values = vec![];
        values.extend(Comm::get_values(trace));
        values.extend(Comm::get_values(quotient_chunks));
        if let Some(random) = random {
            values.extend(Comm::get_values(random));
        }
        values
    }

    fn set_wires(
        &self,
        circuit: &mut CircuitBuilder<F, D>,
        input: Self::Input,
    ) -> Result<(), crate::circuit_builder::CircuitError> {
        let Commitments {
            trace,
            quotient_chunks,
            random,
        } = input;

        let CommitmentWires {
            trace_wires,
            quotient_chunks_wires,
            random_commit,
            _phantom,
        } = self;

        Comm::set_wires(trace_wires, circuit, trace)?;
        Comm::set_wires(quotient_chunks_wires, circuit, quotient_chunks)?;
        if let Some(random_commit) = random_commit {
            let r = random.ok_or(CircuitError::RandomizationError)?;
            Comm::set_wires(random_commit, circuit, r)?;
        }

        Ok(())
    }

    fn num_challenges(&self) -> usize {
        0
    }
}

impl<SC: StarkGenericConfig, const D: usize> Recursive<Val<SC>, D> for OpenedValuesWires<SC, D> {
    type Input = OpenedValues<SC::Challenge>;

    fn new(
        circuit: &mut CircuitBuilder<Val<SC>, D>,
        lens: &mut impl Iterator<Item = usize>,
    ) -> Self {
        let trace_local_len = lens.next().unwrap();
        let mut trace_local_wires = Vec::with_capacity(trace_local_len);
        for _ in 0..trace_local_len {
            trace_local_wires.push(circuit.new_extension_wires());
        }
        let trace_next_len = lens.next().unwrap();
        let mut trace_next_wires = Vec::with_capacity(trace_next_len);
        for _ in 0..trace_next_len {
            trace_next_wires.push(circuit.new_extension_wires());
        }
        let quotient_chunks_len = lens.next().unwrap();
        let quotient_chunks_cols_len = lens.next().unwrap();
        let mut quotient_chunks_wires = Vec::with_capacity(quotient_chunks_len);
        for _ in 0..quotient_chunks_len {
            let mut quotient_col = Vec::with_capacity(quotient_chunks_cols_len);
            for _ in 0..quotient_chunks_cols_len {
                quotient_col.push(circuit.new_extension_wires());
            }
            quotient_chunks_wires.push(quotient_col);
        }
        let random_len = lens.next().unwrap();
        let random_wires = if random_len > 0 {
            let mut r = Vec::with_capacity(random_len);
            for _ in 0..random_len {
                r.push(circuit.new_extension_wires());
            }
            Some(r)
        } else {
            None
        };

        Self {
            trace_local_wires,
            trace_next_wires,
            quotient_chunks_wires,
            random_wires,
            _phantom: PhantomData,
        }
    }

    fn get_values(input: Self::Input) -> Vec<Val<SC>> {
        let OpenedValues {
            trace_local,
            trace_next,
            quotient_chunks,
            random,
        } = input;

        let mut values = vec![];
        values.extend(
            trace_local
                .iter()
                .flat_map(|t| t.as_basis_coefficients_slice()),
        );
        values.extend(
            trace_next
                .iter()
                .flat_map(|t| t.as_basis_coefficients_slice()),
        );
        for chunk in quotient_chunks {
            values.extend(chunk.iter().flat_map(|t| t.as_basis_coefficients_slice()));
        }
        if let Some(random) = random {
            values.extend(random.iter().flat_map(|t| t.as_basis_coefficients_slice()));
        }

        values
    }

    fn set_wires(
        &self,
        circuit: &mut CircuitBuilder<Val<SC>, D>,
        input: Self::Input,
    ) -> Result<(), CircuitError> {
        let OpenedValues {
            trace_local,
            trace_next,
            quotient_chunks,
            random,
        } = input;

        let OpenedValuesWires {
            trace_local_wires,
            trace_next_wires,
            quotient_chunks_wires,
            random_wires,
            ..
        } = self;

        for (wire, value) in trace_local_wires.iter().zip(trace_local.iter()) {
            circuit.set_extension_wires(
                *wire,
                value.as_basis_coefficients_slice().try_into().unwrap(),
            )?;
        }
        for (wire, value) in trace_next_wires.iter().zip(trace_next.iter()) {
            circuit.set_extension_wires(
                *wire,
                value.as_basis_coefficients_slice().try_into().unwrap(),
            )?;
        }
        for (w_chunk, chunk) in quotient_chunks_wires.iter().zip(quotient_chunks.iter()) {
            for (wire, value) in w_chunk.iter().zip(chunk) {
                circuit.set_extension_wires(
                    *wire,
                    value.as_basis_coefficients_slice().try_into().unwrap(),
                )?;
            }
        }
        if let Some(r_wires) = random_wires {
            for (r_wire, r) in r_wires
                .iter()
                .zip_eq(random.ok_or_else(|| CircuitError::RandomizationError)?)
            {
                circuit.set_extension_wires(
                    *r_wire,
                    r.as_basis_coefficients_slice().try_into().unwrap(),
                )?;
            }
        }

        Ok(())
    }

    fn num_challenges(&self) -> usize {
        0
    }
}
