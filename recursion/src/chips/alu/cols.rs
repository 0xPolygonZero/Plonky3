use std::borrow::{Borrow, BorrowMut};
use std::marker::PhantomData;

use p3_field::Field;

use crate::chips::ext_alu::binomial_extension::BinomialExtension;

pub struct AddEvent<F, const R: usize = 1>(pub FieldOpEvent<usize, F, R>);
pub struct SubEvent<F, const R: usize = 1>(pub FieldOpEvent<usize, F, R>);
pub struct MulEvent<F, const R: usize = 1>(pub FieldOpEvent<usize, F, R>);

pub struct ExtAddEvent<F, const D: usize, const R: usize = 1>(
    pub FieldOpEvent<[usize; D], BinomialExtension<F, D>, R>,
);
pub struct ExtSubEvent<F, const D: usize, const R: usize = 1>(
    pub FieldOpEvent<[usize; D], BinomialExtension<F, D>, R>,
);
pub struct ExtMulEvent<F, const D: usize, const R: usize = 1>(
    pub FieldOpEvent<[usize; D], BinomialExtension<F, D>, R>,
);

#[derive(Debug)]
/// Represents an event in the field operation trace.
pub struct FieldOpEvent<W, T, const R: usize = 1> {
    pub left_addr: [W; R],
    pub left_val: [T; R],
    pub right_addr: [W; R],
    pub right_val: [T; R],
    pub res_addr: [W; R],
    pub res_val: [T; R],
}

#[repr(C)]
/// Represents the columns in the ALU trace.
/// R counts how many `a * b = c` operations to do per row in the AIR
pub struct AluCols<F, const R: usize = 1> {
    pub left_addr: [F; R],
    pub left_val: [F; R],
    pub right_addr: [F; R],
    pub right_val: [F; R],
    pub res_addr: [F; R],
    pub res_val: [F; R],
    _phantom_data: PhantomData<F>,
}

impl<F, const R: usize> AluCols<F, R> {
    pub const TRACE_WIDTH: usize = 6 * R;
}

impl<F, const R: usize> Borrow<AluCols<F, R>> for [F] {
    fn borrow(&self) -> &AluCols<F, R> {
        debug_assert_eq!(self.len(), AluCols::<F, R>::TRACE_WIDTH);
        let (prefix, shorts, _suffix) = unsafe { self.align_to::<AluCols<F, R>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<F: Field, const R: usize> BorrowMut<AluCols<F, R>> for [F] {
    fn borrow_mut(&mut self) -> &mut AluCols<F, R> {
        debug_assert_eq!(self.len(), AluCols::<F, R>::TRACE_WIDTH);
        let (prefix, shorts, _suffix) = unsafe { self.align_to_mut::<AluCols<F, R>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}
