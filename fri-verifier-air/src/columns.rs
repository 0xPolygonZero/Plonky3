use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

/// Column layout for FRI commit phase verification.
///
/// Each row represents one FRI folding step with all intermediate arithmetic values
/// needed to verify the folding computation: f(β) = eval₀ + (β - x₀) × (eval₁ - eval₀) × (x₁ - x₀)⁻¹
#[repr(C)]
pub struct CommitPhaseCols<T, const DIGEST_ELEMS: usize = 8> {
    pub query_index: T,
    pub phase_index: T,

    pub beta: [T; 4],
    pub eval_0: [T; 4],
    pub eval_1: [T; 4], 
    pub roll_in_value: [T; 4],

    pub domain_index: T,
    pub parent_index: T,
    pub sibling_index: T,

    pub x0: T,
    pub x1: T,
    pub mmcs_verified: T,
    pub log_height: T,
    
    // FRI folding intermediate arithmetic
    pub beta_minus_x0: [T; 4],
    pub eval_diff: [T; 4],
    pub x_diff: T,
    pub x_diff_inv: T,
    pub beta_eval_product: [T; 4],
    pub interpolation_term: [T; 4],
    pub folded_eval: [T; 4],
    pub beta_squared: [T; 4],
    pub roll_in_contribution: [T; 4],

    // MMCS verification
    pub fri_commit: [T; DIGEST_ELEMS],
    pub sibling_eval: [T; 4],
}

pub const fn num_commit_phase_cols<const DIGEST_ELEMS: usize>() -> usize {
    size_of::<CommitPhaseCols<u8, DIGEST_ELEMS>>()
}

// Default version for backward compatibility
pub const fn num_commit_phase_cols_default() -> usize {
    num_commit_phase_cols::<8>()
}

impl<T, const DIGEST_ELEMS: usize> Borrow<CommitPhaseCols<T, DIGEST_ELEMS>> for [T] {
    fn borrow(&self) -> &CommitPhaseCols<T, DIGEST_ELEMS> {
        debug_assert_eq!(self.len(), num_commit_phase_cols::<DIGEST_ELEMS>());
        let (prefix, shorts, suffix) = unsafe { self.align_to::<CommitPhaseCols<T, DIGEST_ELEMS>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T, const DIGEST_ELEMS: usize> BorrowMut<CommitPhaseCols<T, DIGEST_ELEMS>> for [T] {
    fn borrow_mut(&mut self) -> &mut CommitPhaseCols<T, DIGEST_ELEMS> {
        debug_assert_eq!(self.len(), num_commit_phase_cols::<DIGEST_ELEMS>());
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<CommitPhaseCols<T, DIGEST_ELEMS>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}
