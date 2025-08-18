use alloc::vec;
use alloc::vec::Vec;
use itertools::Itertools;
use p3_field::{ExtensionField, Field};

use crate::{Entry, SymbolicExpression};

#[derive(Debug, Clone, Copy)]
pub enum Gate {
    Add,
    Mul,
    Sub,
    Fri,
    MerkleTree,
}

#[derive(Copy, Clone, Debug)]
pub struct Wire<F: Field> {
    pub value: Option<F>,
    pub index: usize, // Index in the witness table
}

pub trait CircuitBuilder<F: Field> {
    fn add_gate(&mut self, g: Gate, a: Wire<F>, b: Wire<F>) -> Wire<F>;

    fn new_wire(&mut self) -> Wire<F>;

    fn add_constant(&mut self, value: F) -> Wire<F>;
}

pub fn symbolic_to_circuit<F: Field, EF: ExtensionField<F>, C: CircuitBuilder<EF>>(
    is_first_row: Wire<EF>,
    is_last_row: Wire<EF>,
    is_transition: Wire<EF>,
    challenges: &[Wire<EF>],
    public_values: &[Wire<EF>],
    local_prep_values: &[Wire<EF>],
    next_prep_values: &[Wire<EF>],
    local_values: &[Wire<EF>],
    next_values: &[Wire<EF>],
    symbolic: &SymbolicExpression<F>,
    circuit: &mut C,
) -> Wire<EF> {
    match symbolic {
        SymbolicExpression::Constant(c) => circuit.add_constant(EF::from(c.clone())),
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
            Entry::Public => public_values[v.index].clone(),
            Entry::Challenge => challenges[v.index].clone(),
            _ => unimplemented!(),
        },
        SymbolicExpression::Add { x, y, .. } => {
            let x_wire = symbolic_to_circuit(
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
            let y_wire = symbolic_to_circuit(
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
            circuit.add_gate(Gate::Add, x_wire, y_wire)
        }
        SymbolicExpression::Mul { x, y, .. } => {
            let x_wire = symbolic_to_circuit(
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
            let y_wire = symbolic_to_circuit(
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
            circuit.add_gate(Gate::Mul, x_wire, y_wire)
        }
        SymbolicExpression::Sub { x, y, .. } => {
            let x_wire = symbolic_to_circuit(
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
            let y_wire = symbolic_to_circuit(
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

            circuit.add_gate(Gate::Sub, x_wire, y_wire)
        }
        SymbolicExpression::Neg { x, .. } => {
            let x_wire = symbolic_to_circuit(
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
            let zero = circuit.add_constant(EF::ZERO);

            circuit.add_gate(Gate::Sub, zero, x_wire)
        }
        SymbolicExpression::IsFirstRow => is_first_row,
        SymbolicExpression::IsLastRow => is_last_row,
        SymbolicExpression::IsTransition => is_transition,
    }
}

pub struct DummyCircuit<F: Field> {
    pub wires: Vec<Wire<F>>,
    pub events: Vec<(Gate, Wire<F>, Wire<F>, Wire<F>)>,
    forest: Vec<Vec<(usize, usize)>>,
}

impl<F: Field> CircuitBuilder<F> for DummyCircuit<F> {
    fn add_gate(&mut self, g: Gate, a: Wire<F>, b: Wire<F>) -> Wire<F> {
        let output = self.new_wire();
        self.events.push((g, a.clone(), b.clone(), output.clone()));

        let output_index = self.wires.len() - 1;
        let a_index = a.index;
        let b_index = b.index;

        assert!(self.forest.len() == output.index + 1);

        // Update the forest
        self.forest[a_index].push((self.events.len() - 1, 0));
        self.forest[b_index].push((self.events.len() - 1, 1));
        self.forest[output_index].push((self.events.len() - 1, 2));
        output
    }

    fn new_wire(&mut self) -> Wire<F> {
        let index = self.wires.len();
        let wire = Wire { value: None, index };
        self.wires.push(wire.clone());

        self.forest.push(vec![]);
        wire
    }

    fn add_constant(&mut self, value: F) -> Wire<F> {
        let mut wire = self.new_wire();
        wire.value = Some(value);
        wire
    }
}

impl<F: Field> DummyCircuit<F> {
    pub fn new() -> Self {
        DummyCircuit {
            wires: Vec::new(),
            events: Vec::new(),
            forest: Vec::new(),
        }
    }

    fn set_wire(&mut self, value: F, wire: &mut Wire<F>) {
        wire.value = Some(value);

        for (row, col) in &self.forest[wire.index] {
            match col {
                0 => self.events[*row].1.value = Some(value),
                1 => self.events[*row].2.value = Some(value),
                2 => self.events[*row].3.value = Some(value),
                _ => unreachable!(),
            }
        }
    }

    fn set_wire_by_index(&mut self, value: F, wire_idx: usize) {
        self.wires[wire_idx].value = Some(value);

        for (row, col) in &self.forest[wire_idx].clone() {
            match col {
                0 => self.events[*row].1.value = Some(value),
                1 => self.events[*row].2.value = Some(value),
                2 => self.events[*row].3.value = Some(value),
                _ => unreachable!(),
            }
        }
    }

    pub fn generate(
        &mut self,
        sels: &[F],
        public_inputs: &[F],
        local_values: &[F],
        next_values: &[F],
        sels_circuit: &mut [Wire<F>],
        public_inputs_circuit: &mut [Wire<F>],
        local_values_circuit: &mut [Wire<F>],
        next_values_circuit: &mut [Wire<F>],
    ) {
        for (s_c, s) in sels_circuit.iter_mut().zip_eq(sels.iter()) {
            self.set_wire(*s, s_c);
        }

        for (p_i_c, p_i) in public_inputs_circuit
            .iter_mut()
            .zip_eq(public_inputs.iter())
        {
            self.set_wire(*p_i, p_i_c);
        }

        for (l_v_c, l_v) in local_values_circuit.iter_mut().zip_eq(local_values.iter()) {
            self.set_wire(*l_v, l_v_c);
        }

        for (n_v_c, n_v) in next_values_circuit.iter_mut().zip_eq(next_values.iter()) {
            self.set_wire(*n_v, n_v_c);
        }

        for i in 0..self.events.len() {
            let event = &self.events[i];
            let (out, out_wire_idx) = match event.0 {
                Gate::Add => {
                    let out = event.1.value.unwrap() + event.2.value.unwrap();
                    (out, event.3.index)
                }
                Gate::Mul => {
                    let out = event.1.value.unwrap() * event.2.value.unwrap();
                    (out, event.3.index)
                }
                Gate::Sub => {
                    let out = event.1.value.unwrap() - event.2.value.unwrap();
                    (out, event.3.index)
                }
                Gate::Fri | Gate::MerkleTree => {
                    // These gates do not produce a value in this dummy circuit
                    unimplemented!()
                }
            };

            // Now we can safely call set_wire_by_index
            self.set_wire_by_index(out, out_wire_idx);
        }
    }
}
