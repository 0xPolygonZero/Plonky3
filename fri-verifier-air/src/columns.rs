use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

/// Columns for the CommitPhase table of FRI verifier AIR.
///
/// Each row represents one FRI folding step in the commit phase.
/// Maps to the FRI folding logic in verifier.rs.
#[repr(C)]
pub struct CommitPhaseCols<T> {
    // Cross-table lookup to QueryRound (base field)
    /// Links to QueryRound.query_index
    pub query_index: T,

    // Phase identification (base field)
    /// Phase index (0 to commit_phase_commits.len()-1)
    pub phase_index: T,

    // FRI challenges and evaluations (extension field elements stored as [T; 4])
    /// Beta challenge from betas vector (4 components)
    pub beta: [T; 4],

    /// eval_0 = folded_eval from previous phase (4 components)
    pub eval_0: [T; 4],

    /// eval_1 = opening.sibling_value (4 components)
    pub eval_1: [T; 4],

    /// roll_in_value = beta^2 * reduced_opening (4 components)
    pub roll_in_value: [T; 4],

    // Index management during folding (base field)
    /// Current index before folding  
    pub domain_index: T,
    /// domain_index >> 1 (after folding)
    pub parent_index: T,
    /// domain_index ^ 1
    pub sibling_index: T,

    // Subgroup points for interpolation (base field)
    /// subgroup_start  
    pub x0: T,
    /// -subgroup_start
    pub x1: T,

    // Verification flags (base field)
    /// MMCS verification result from verifier.rs
    pub mmcs_verified: T,
    /// Current domain height: log_max_height - phase_index - 1
    pub log_height: T,
    // NOTE: Extension field arithmetic operations (beta - x0, eval_1 - eval_0, etc.)
    // will be handled by dedicated extension field operation tables.
    // This table focuses on the high-level FRI folding logic and uses
    // cross-table lookups to verify the extension field computations.
}

pub const fn num_commit_phase_cols() -> usize {
    size_of::<CommitPhaseCols<u8>>()
}

impl<T> Borrow<CommitPhaseCols<T>> for [T] {
    fn borrow(&self) -> &CommitPhaseCols<T> {
        debug_assert_eq!(self.len(), num_commit_phase_cols());
        let (prefix, shorts, suffix) = unsafe { self.align_to::<CommitPhaseCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T> BorrowMut<CommitPhaseCols<T>> for [T] {
    fn borrow_mut(&mut self) -> &mut CommitPhaseCols<T> {
        debug_assert_eq!(self.len(), num_commit_phase_cols());
        let (prefix, shorts, suffix) = unsafe { self.align_to_mut::<CommitPhaseCols<T>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}
