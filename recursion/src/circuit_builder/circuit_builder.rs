use std::array;

use itertools::Itertools;
use p3_commit::Pcs;
use p3_field::extension::BinomiallyExtendable;
use p3_field::{BasedVectorSpace, ExtensionField, Field};
use p3_uni_stark::{
    Commitments, Entry, OpenedValues, Proof, StarkGenericConfig, SymbolicExpression, Val,
};

use crate::chips::witness::air::RomEvent;
use crate::circuit_builder::gates::arith_gates::{
    AddExtensionGate, MulExtensionGate, SubExtensionGate,
};
use crate::circuit_builder::gates::event::AllEvents;
use crate::circuit_builder::gates::gate::Gate;
use crate::verifier::circuit_verifier::{CommitmentWires, OpenedValuesWires, ProofWires};
use crate::verifier::recursive_traits::{ForRecursiveVersion, RecursiveVersion};

pub type WireId = usize;

#[derive(Default)]
pub struct CircuitBuilder<F: Field, const D: usize> {
    wires: Vec<Option<F>>,
    gate_instances: Vec<Box<dyn Gate<F, D>>>,
}

pub type ExtensionWireId<const D: usize> = [WireId; D];

impl<F: Field, const D: usize> CircuitBuilder<F, D> {
    pub fn new() -> Self {
        CircuitBuilder {
            wires: Vec::new(),
            gate_instances: Vec::new(),
            ..Default::default()
        }
    }

    pub fn new_wire(&mut self) -> WireId {
        self.wires.push(None);
        self.wires.len() - 1
    }

    pub fn new_extension_wires(&mut self) -> ExtensionWireId<D> {
        array::from_fn(|_| self.new_wire())
    }

    pub fn add_gate(&mut self, gate: Box<dyn Gate<F, D>>) {
        self.gate_instances.push(gate);
    }

    pub fn add_constant(&mut self, value: F) -> WireId {
        let wire_id = self.new_wire();
        self.set_wire_value(wire_id, value).unwrap();
        wire_id
    }

    pub fn add_extension_constant<EF: ExtensionField<F>>(
        &mut self,
        value: EF,
    ) -> ExtensionWireId<D> {
        assert_eq!(
            EF::DIMENSION,
            D,
            "The dimension of the extension field must match the number of wires in a ChallengerWire"
        );
        let wire_id = array::from_fn(|_| self.new_wire());
        let base_values = value.as_basis_coefficients_slice();
        self.set_wire_values(&wire_id, &base_values).unwrap();
        wire_id
    }

    pub fn wires(&mut self) -> &Vec<Option<F>> {
        &mut self.wires
    }

    pub fn get_wire_value(&self, id: WireId) -> Result<Option<F>, CircuitError> {
        if id >= self.wires.len() {
            Err(CircuitError::InvalidWireId)
        } else {
            Ok(self.wires[id])
        }
    }

    pub fn set_wire_value(&mut self, id: WireId, value: F) -> Result<(), CircuitError> {
        let prev = self.get_wire_value(id)?;
        if let Some(val) = prev {
            if val != value {
                return Err(CircuitError::WireSetTwice);
            }
        } else {
            self.wires[id] = Some(value);
        }
        Ok(())
    }

    pub fn set_wire_values(&mut self, ids: &[WireId], value: &[F]) -> Result<(), CircuitError> {
        for (id, val) in ids.iter().zip(value.iter()) {
            self.set_wire_value(*id, *val)?;
        }

        Ok(())
    }

    pub fn set_extension_wires(
        &mut self,
        ids: ExtensionWireId<D>,
        value: &[F; D],
    ) -> Result<(), CircuitError> {
        self.set_wire_values(&ids, value)
    }

    pub fn set_extension_wires_slice(
        &mut self,
        ids: &[ExtensionWireId<D>],
        value: &[[F; D]],
    ) -> Result<(), CircuitError> {
        for (id, val) in ids.iter().zip(value.iter()) {
            self.set_extension_wires(*id, val)?;
        }
        Ok(())
    }

    pub fn set_public_inputs(&mut self, ids: &[WireId], values: &[F]) -> Result<(), CircuitError> {
        self.set_wire_values(ids, values)
    }

    pub fn get_instances(&self) -> &Vec<Box<dyn Gate<F, D>>> {
        &self.gate_instances
    }

    pub fn generate(&mut self) -> Result<AllEvents<F, D>, CircuitError> {
        let mut gate_instances = core::mem::take(&mut self.gate_instances);
        let mut all_events = AllEvents {
            add_events: Vec::new(),
            sub_events: Vec::new(),
            mul_events: Vec::new(),
            ext_add_events: Vec::new(),
            ext_sub_events: Vec::new(),
            ext_mul_events: Vec::new(),
            witness_events: Vec::new(),
        };

        // Generaete events for all gates
        for gate in gate_instances.iter_mut() {
            gate.generate(self, &mut all_events)?;
        }

        // Generate the witeness events from the wire values
        all_events.witness_events = self
            .wires
            .iter()
            .enumerate()
            .filter_map(|(i, v)| v.map(|val| RomEvent(i, val)))
            .collect(); // TODO: avoid the intermediate allocation?

        Ok(all_events)
    }
}

/// Function which, given `CommitmentWires` and `Commitments`, sets the wires to the associated values.
pub fn set_commitment_wires<SC: StarkGenericConfig, Comm: RecursiveVersion, const D: usize>(
    circuit: &mut CircuitBuilder<Val<SC>, D>,
    comm_wires: &CommitmentWires<Comm>,
    comm: &Commitments<<SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment>,
) -> Result<(), CircuitError>
where
    <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment: ForRecursiveVersion<Val<SC>>,
{
    let CommitmentWires {
        trace_wires: trace_wires_comm,
        quotient_chunks_wires: quotient_chunks_wires_comm,
        random_commit,
        ..
    } = comm_wires;

    circuit.set_wire_values(&trace_wires_comm.get_wires(), &comm.trace.get_values())?;
    circuit.set_wire_values(
        &quotient_chunks_wires_comm.get_wires(),
        &comm.quotient_chunks.get_values(),
    )?;
    if let Some(r_commit) = random_commit {
        circuit.set_wire_values(
            &r_commit.get_wires(),
            &comm
                .random
                .as_ref()
                .ok_or(CircuitError::RandomizationError)?
                .get_values(),
        )?;
    }
    Ok(())
}

/// Function which, given `OpenedValuesWires` and `OpenedValues`, sets the wires to the aassociated values.
pub fn set_opened_wires<SC: StarkGenericConfig, Comm: RecursiveVersion, const D: usize>(
    circuit: &mut CircuitBuilder<Val<SC>, D>,
    opened_wires: &OpenedValuesWires<D>,
    opened_values: OpenedValues<SC::Challenge>,
) -> Result<(), CircuitError>
where
    <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment: ForRecursiveVersion<Val<SC>>,
{
    let OpenedValuesWires {
        trace_local_wires,
        trace_next_wires,
        quotient_chunks_wires,
        random,
    } = opened_wires;

    let trace_local_ext = opened_values
        .trace_local
        .iter()
        .map(|l| l.as_basis_coefficients_slice().try_into().unwrap())
        .collect::<Vec<_>>();
    circuit.set_extension_wires_slice(trace_local_wires, &trace_local_ext)?;
    let trace_next_ext = opened_values
        .trace_next
        .iter()
        .map(|l| l.as_basis_coefficients_slice().try_into().unwrap())
        .collect::<Vec<_>>();
    circuit.set_extension_wires_slice(trace_next_wires, &trace_next_ext)?;
    for (chunk_wires, chunk) in quotient_chunks_wires
        .iter()
        .zip_eq(opened_values.quotient_chunks.iter())
    {
        let chunk_ext = chunk
            .iter()
            .map(|l| l.as_basis_coefficients_slice().try_into().unwrap())
            .collect::<Vec<_>>();
        circuit.set_extension_wires_slice(chunk_wires, &chunk_ext)?;
    }

    if let Some(r) = opened_values.random {
        let r_ext = r
            .iter()
            .map(|l| l.as_basis_coefficients_slice().try_into().unwrap())
            .collect::<Vec<_>>();
        let random_wires = random.as_ref().ok_or(CircuitError::RandomizationError)?;
        circuit.set_extension_wires_slice(random_wires, &r_ext)?;
    }
    Ok(())
}

/// Given a proof and proof wires, this function sets the proof wires with the corresponding values.
pub fn set_proof_wires<
    SC: StarkGenericConfig,
    Comm: RecursiveVersion,
    InputProof: RecursiveVersion,
    const D: usize,
>(
    circuit: &mut CircuitBuilder<Val<SC>, D>,
    proof_wires: &ProofWires<D, Comm, InputProof>,
    proof: Proof<SC>,
) -> Result<(), CircuitError>
where
    <SC::Pcs as Pcs<SC::Challenge, SC::Challenger>>::Commitment: ForRecursiveVersion<Val<SC>>,
    <<SC as StarkGenericConfig>::Pcs as Pcs<
        <SC as StarkGenericConfig>::Challenge,
        <SC as StarkGenericConfig>::Challenger,
    >>::Proof: ForRecursiveVersion<Val<SC>>,
{
    let ProofWires {
        commitments_wires,
        opened_values_wires,
        opening_proof: opening_proof_wires,
        ..
    } = proof_wires;

    let Proof {
        commitments,
        opened_values,
        opening_proof,
        degree_bits: _,
    } = proof;

    // Set commitment wires.
    set_commitment_wires::<SC, Comm, D>(circuit, commitments_wires, &commitments)?;

    // Set opened values.
    set_opened_wires::<SC, Comm, D>(circuit, opened_values_wires, opened_values)?;

    let opening_proof_wires = opening_proof_wires.get_wires();
    let opening_proof_values = opening_proof.get_values();
    circuit.set_wire_values(&opening_proof_wires, &opening_proof_values)?;
    Ok(())
}

#[derive(Debug)]
pub enum CircuitError {
    RandomizationError,
    InvalidProofShape,
    InvalidWireId,
    InputNotSet,
    WireSetTwice,
}

pub fn symbolic_to_circuit<F: Field, EF, const D: usize>(
    is_first_row: ExtensionWireId<D>,
    is_last_row: ExtensionWireId<D>,
    is_transition: ExtensionWireId<D>,
    challenges: &[ExtensionWireId<D>],
    public_values: &[WireId],
    local_prep_values: &[ExtensionWireId<D>],
    next_prep_values: &[ExtensionWireId<D>],
    local_values: &[ExtensionWireId<D>],
    next_values: &[ExtensionWireId<D>],
    symbolic: &SymbolicExpression<EF>,
    circuit: &mut CircuitBuilder<F, D>,
) -> ExtensionWireId<D>
where
    F: BinomiallyExtendable<D>,
    EF: ExtensionField<F>,
{
    assert_eq!(D, EF::DIMENSION);
    match symbolic {
        SymbolicExpression::Constant(c) => circuit.add_extension_constant(c.clone()),
        SymbolicExpression::Variable(v) => match v.entry {
            Entry::Preprocessed { offset } => {
                if offset == 0 {
                    local_prep_values[v.index].clone()
                } else if offset == 1 {
                    next_prep_values[v.index].clone()
                } else {
                    panic!("Cannot have expressions involving more than two rows.")
                }
            }
            Entry::Main { offset } => {
                if offset == 0 {
                    local_values[v.index].clone()
                } else if offset == 1 {
                    next_values[v.index].clone()
                } else {
                    panic!("Cannot have expressions involving more than two rows.")
                }
            }
            Entry::Public => {
                let mut res = circuit.add_extension_constant(EF::ZERO);
                res[0] = public_values[v.index].clone();

                res
            }
            Entry::Challenge => challenges[v.index].clone(),
            _ => unimplemented!(),
        },
        SymbolicExpression::Add { x, y, .. } => {
            let x_wire = symbolic_to_circuit::<F, EF, D>(
                is_first_row.clone(),
                is_last_row.clone(),
                is_transition.clone(),
                challenges,
                public_values,
                local_prep_values,
                next_prep_values,
                local_values,
                next_values,
                x,
                circuit,
            );
            let y_wire = symbolic_to_circuit::<F, EF, D>(
                is_first_row,
                is_last_row,
                is_transition,
                challenges,
                public_values,
                local_prep_values,
                next_prep_values,
                local_values,
                next_values,
                y,
                circuit,
            );

            let out_wire = circuit.new_extension_wires();

            AddExtensionGate::add_to_circuit(circuit, x_wire, y_wire, out_wire);

            out_wire
        }
        SymbolicExpression::Mul { x, y, .. } => {
            let x_wire = symbolic_to_circuit::<F, EF, D>(
                is_first_row.clone(),
                is_last_row.clone(),
                is_transition.clone(),
                challenges,
                public_values,
                local_prep_values,
                next_prep_values,
                local_values,
                next_values,
                x,
                circuit,
            );
            let y_wire = symbolic_to_circuit::<F, EF, D>(
                is_first_row,
                is_last_row,
                is_transition,
                challenges,
                public_values,
                local_prep_values,
                next_prep_values,
                local_values,
                next_values,
                y,
                circuit,
            );

            let out_wire = circuit.new_extension_wires();

            MulExtensionGate::add_to_circuit(circuit, x_wire, y_wire, out_wire);
            out_wire
        }
        SymbolicExpression::Sub { x, y, .. } => {
            let x_wire = symbolic_to_circuit::<F, EF, D>(
                is_first_row.clone(),
                is_last_row.clone(),
                is_transition.clone(),
                challenges,
                public_values,
                local_prep_values,
                next_prep_values,
                local_values,
                next_values,
                x,
                circuit,
            );
            let y_wire = symbolic_to_circuit::<F, EF, D>(
                is_first_row,
                is_last_row,
                is_transition,
                challenges,
                public_values,
                local_prep_values,
                next_prep_values,
                local_values,
                next_values,
                y,
                circuit,
            );

            let out_wire = circuit.new_extension_wires();

            SubExtensionGate::add_to_circuit(circuit, x_wire, y_wire, out_wire);

            out_wire
        }
        SymbolicExpression::Neg { x, .. } => {
            let x_wire = symbolic_to_circuit::<F, EF, D>(
                is_first_row,
                is_last_row,
                is_transition,
                challenges,
                public_values,
                local_prep_values,
                next_prep_values,
                local_values,
                next_values,
                x,
                circuit,
            );
            let zero = circuit.add_extension_constant(EF::ZERO);

            let out_wire = circuit.new_extension_wires();

            SubExtensionGate::add_to_circuit(circuit, zero, x_wire, out_wire);

            out_wire
        }
        SymbolicExpression::IsFirstRow => is_first_row,
        SymbolicExpression::IsLastRow => is_last_row,
        SymbolicExpression::IsTransition => is_transition,
    }
}
