use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

/// Columns for the CommitPhase table of FRI verifier AIR.
///
/// Each row represents one FRI folding step in the commit phase.
/// Maps to the folding loop at fri.rs:727-908 in simulate_fri_verification_immediate_values.
#[repr(C)]
pub struct CommitPhaseCols<T> {
    // Cross-table lookup to QueryRound
    /// Links to QueryRound.query_index
    pub query_id: T,

    // Phase identification (fri.rs:727)
    /// Phase index (0 to commit_phase_commits.len()-1)
    pub phase_idx: T,

    // Challenge for this phase (fri.rs:727)
    /// Beta challenge from betas vector
    pub beta: T,

    // Index management during folding (fri.rs:766)
    /// Current index before folding  
    pub domain_index: T,
    /// domain_index >> 1 (after folding)
    pub parent_index: T,
    /// domain_index ^ 1
    pub sibling_index: T,

    // The two evaluations being folded (fri.rs:754-755)
    /// evals[0] = folded_eval from previous phase
    pub eval_0: T,
    /// evals[1] = opening.sibling_value
    pub eval_1: T,

    // Subgroup points for interpolation (fri.rs:847-848)
    /// subgroup_start  
    pub x0: T,
    /// -subgroup_start
    pub x1: T,

    // FRI folding computation (fri.rs:857-874)
    /// beta - x0
    pub beta_minus_x0: T,
    /// eval_1 - eval_0  
    pub eval_1_minus_eval_0: T,
    /// x1 - x0
    pub x1_minus_x0: T,
    /// (x1 - x0)^(-1)
    pub inverse: T,
    /// (beta - x0) * (eval_1 - eval_0) * (x1 - x0)^(-1)
    pub intermediate: T,
    /// eval_0 + intermediate (before roll-in)
    pub folded_result_pre: T,

    // Roll-in from reduced openings (fri.rs:883-903)
    /// beta^2 * reduced_opening (0 if no roll-in)
    pub roll_in_value: T,
    /// Final folded evaluation after roll-in
    pub folded_result: T,

    // Verification flags
    /// MMCS verification result (fri.rs:791-807)
    pub mmcs_verified: T,
    /// Current domain height: log_max_height - phase_idx - 1
    pub log_height: T,
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
