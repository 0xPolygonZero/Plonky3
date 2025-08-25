use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;

use crate::air::alu::cols::FieldOpEvent;
use crate::air::witness::air::RomAir;
use crate::circuit_builder::gates::event::AllEvents;

pub trait Table<F: Field> {
    fn generate_trace(&self, all_events: &AllEvents<F>) -> RowMajorMatrix<F>;
}

pub struct AddTable;
pub struct SubTable;
pub struct MulTable;

pub struct WitnessTable;

// TODO: Maybe do this with a macro?
impl<F: Field> Table<F> for AddTable {
    fn generate_trace(&self, all_events: &AllEvents<F>) -> RowMajorMatrix<F> {
        FieldOpEvent::generate_trace(
            all_events.add_events.iter().map(|x| &x.0),
            all_events.add_events.len(),
        )
    }
}

impl<F: Field> Table<F> for SubTable {
    fn generate_trace(&self, all_events: &AllEvents<F>) -> RowMajorMatrix<F> {
        FieldOpEvent::generate_trace(
            all_events.sub_events.iter().map(|x| &x.0),
            all_events.sub_events.len(),
        )
    }
}

impl<F: Field> Table<F> for MulTable {
    fn generate_trace(&self, all_events: &AllEvents<F>) -> RowMajorMatrix<F> {
        FieldOpEvent::generate_trace(
            all_events.mul_events.iter().map(|x| &x.0),
            all_events.mul_events.len(),
        )
    }
}

impl<F: Field> Table<F> for WitnessTable {
    fn generate_trace(&self, all_events: &AllEvents<F>) -> RowMajorMatrix<F> {
        RomAir::build_trace(
            all_events.witness_events.iter(),
            all_events.witness_events.len(),
        )
    }
}
