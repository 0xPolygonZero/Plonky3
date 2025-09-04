use p3_field::Field;

use crate::chips::alu::cols::{ExtAddEvent, ExtMulEvent, ExtSubEvent};
use crate::chips::merkle_path::air::MerklePathEvent;
use crate::chips::witness::air::RomEvent;
use crate::chips::{AddEvent, MulEvent, SubEvent};

#[derive(Default, Debug)]
pub struct AllEvents<F: Field, const D: usize, const DIGEST_ELEMS: usize> {
    pub add_events: Vec<AddEvent<F>>,
    pub sub_events: Vec<SubEvent<F>>,
    pub mul_events: Vec<MulEvent<F>>,

    pub ext_add_events: Vec<ExtAddEvent<F, D>>,
    pub ext_sub_events: Vec<ExtSubEvent<F, D>>,
    pub ext_mul_events: Vec<ExtMulEvent<F, D>>,

    pub merkle_path_events: Vec<MerklePathEvent<F, DIGEST_ELEMS>>,

    pub witness_events: Vec<RomEvent<F>>,
}
