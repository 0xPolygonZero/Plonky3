use std::char::MAX;

use itertools::Itertools;
use p3_field::{Field, PrimeField, PrimeField32};
use p3_symmetric::FieldCompression;

use crate::chips::merkle_path::air::{MerklePathEvent, SiblingWithExtraOpening, number_of_nodes};
use crate::circuit_builder::gates::event::AllEvents;
use crate::circuit_builder::gates::gate::Gate;
use crate::circuit_builder::{CircuitBuilder, CircuitError, WireId};

pub struct MerklePathGate<F, C, const DIGEST_ELEMS: usize, const MAX_TREE_HEIGHT: usize>
where
    F: PrimeField32,
    C: FieldCompression<F, 2, DIGEST_ELEMS> + 'static,
{
    leaf: [WireId; DIGEST_ELEMS],
    index: WireId,
    siblings: Vec<SiblingWithExtraOpening<[WireId; DIGEST_ELEMS]>>,
    root: [WireId; DIGEST_ELEMS],
    compress: C,
    _marker: std::marker::PhantomData<F>,
}

impl<F, C, const DIGEST_ELEMS: usize, const MAX_TREE_HEIGHT: usize>
    MerklePathGate<F, C, DIGEST_ELEMS, MAX_TREE_HEIGHT>
where
    F: PrimeField32,
    C: FieldCompression<F, 2, DIGEST_ELEMS>,
{
    pub fn new(
        leaf: [WireId; DIGEST_ELEMS],
        index: WireId,
        siblings: Vec<SiblingWithExtraOpening<[WireId; DIGEST_ELEMS]>>,
        root: [WireId; DIGEST_ELEMS],
        compress: C,
    ) -> Self {
        debug_assert!(siblings.len() <= MAX_TREE_HEIGHT);

        MerklePathGate {
            leaf,
            index,
            siblings,
            root,
            compress,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn add_to_circuit<const D: usize>(
        builder: &mut CircuitBuilder<F, D, DIGEST_ELEMS>,
        leaf: [WireId; DIGEST_ELEMS],
        index: WireId,
        siblings: Vec<SiblingWithExtraOpening<[WireId; DIGEST_ELEMS]>>,
        root: [WireId; DIGEST_ELEMS],
        compress: C,
    ) {
        let gate = Self::new(leaf, index, siblings, root, compress);
        builder.add_gate(Box::new(gate));
    }
}

impl<F, C, const D: usize, const DIGEST_ELEMS: usize, const MAX_TREE_HEIGHT: usize>
    Gate<F, D, DIGEST_ELEMS> for MerklePathGate<F, C, DIGEST_ELEMS, MAX_TREE_HEIGHT>
where
    F: PrimeField32,
    C: FieldCompression<F, 2, DIGEST_ELEMS>,
{
    fn n_inputs(&self) -> usize {
        self.leaf.len() + DIGEST_ELEMS * self.siblings.iter().map(number_of_nodes).sum::<usize>()
    }

    fn n_outputs(&self) -> usize {
        self.root.len()
    }

    fn generate(
        &self,
        builder: &mut CircuitBuilder<F, D, DIGEST_ELEMS>,
        all_events: &mut AllEvents<F, D, DIGEST_ELEMS>,
    ) -> Result<(), CircuitError> {
        let get_wire_value = |wire| {
            builder
                .get_wire_value(wire)
                .expect("Wire not set")
                .expect("Wire not set") // TOD: Manage errors  properly
        };
        let leaf = self.leaf.map(get_wire_value);
        let index = get_wire_value(self.index);
        let siblings = self
            .siblings
            .iter()
            .map(|(sibling_wires, extra_row_wires)| {
                (
                    sibling_wires.map(get_wire_value),
                    extra_row_wires.map(|extra_row_wires| extra_row_wires.map(get_wire_value)),
                )
            })
            .collect_vec();
        let index = index.as_canonical_u32();
        let index_bits = (0..32).map(|i| (index >> i & 1) == 1);
        let root = siblings.iter().zip(index_bits).fold(
            leaf,
            |state, ((sibling, extra_row), index_bit)| {
                let mut state = if index_bit {
                    self.compress.compress_field([*sibling, state])
                } else {
                    self.compress.compress_field([state, *sibling])
                };
                if let Some(extra_row) = extra_row {
                    state = self.compress.compress_field([state, *extra_row])
                }
                state
            },
        );
        all_events.merkle_path_events.push(MerklePathEvent {
            leaf,
            index,
            siblings,
            root,
        });
        Ok(())
    }
}
