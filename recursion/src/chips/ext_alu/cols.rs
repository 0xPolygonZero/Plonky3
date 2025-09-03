use std::borrow::{Borrow, BorrowMut};

use crate::chips::ext_alu::binomial_extension::BinomialExtension;

#[repr(C)]
/// Represents the columns in the ALU trace.
/// R counts how many `a * b = c` operations to do per row in the AIR
pub struct ExtAluCols<F, const D: usize, const R: usize = 1> {
    pub left_addr: [[F; D]; R],
    pub left_val: [BinomialExtension<F, D>; R],
    pub right_addr: [[F; D]; R],
    pub right_val: [BinomialExtension<F, D>; R],
    pub res_addr: [[F; D]; R],
    pub res_val: [BinomialExtension<F, D>; R],
}

impl<F, const D: usize, const R: usize> ExtAluCols<F, D, R> {
    pub const TRACE_WIDTH: usize = 6 * R * D;
}

impl<F, const R: usize, const D: usize> Borrow<ExtAluCols<F, D, R>> for [F] {
    fn borrow(&self) -> &ExtAluCols<F, D, R> {
        println!("self_len = {:?} self = [", self.len());
        debug_assert_eq!(self.len(), ExtAluCols::<F, D, R>::TRACE_WIDTH);
        let (prefix, shorts, _suffix) = unsafe { self.align_to::<ExtAluCols<F, D, R>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<F, const R: usize, const D: usize> BorrowMut<ExtAluCols<F, D, R>> for [F] {
    fn borrow_mut(&mut self) -> &mut ExtAluCols<F, D, R> {
        debug_assert_eq!(self.len(), ExtAluCols::<F, D, R>::TRACE_WIDTH);
        let (prefix, shorts, _suffix) = unsafe { self.align_to_mut::<ExtAluCols<F, D, R>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}
