use core::borrow::{Borrow, BorrowMut};
use core::mem::size_of;

/// Column layout for FRI commit phase verification.
///
/// Each row represents one FRI folding step with all intermediate arithmetic values
/// needed to verify the folding computation:
///   f(β) = eval₀ + (β - x₀) × (eval₁ - eval₀) × (x₁ - x₀)⁻¹  (+ optional roll-in)
///
/// Notes:
/// - `eval_1` is the sibling evaluation proved against `fri_commit` via MMCS.
/// - `is_roll_in` gates whether a reduced opening rolls in at this height.
/// - `reduced_opening` is fetched via CTL from a ReducedOpenings table keyed by
///    (query_index, log_height) and used only if `is_roll_in = 1`.
/// - Index relations should be enforced by constraints:
///     domain_index = 2 * parent_index + lsb
///     sibling_index = 2 * parent_index + (1 - lsb)
/// - Subgroup points should satisfy x1 = -x0 and x0 ∈ 2-adic subgroup for `log_height+1`.
#[repr(C)]
pub struct CommitPhaseCols<T, const DIGEST_ELEMS: usize = 8> {
    // Cross-table linkages / identifiers (base field)
    /// Links rows across phases for the same query
    pub query_index: T,
    /// Phase index in [0, commit_phase_commits.len())
    pub phase_index: T,

    // FRI challenge & evaluations (extension field split into 4 limbs)
    /// Beta challenge (sampled after observing `fri_commit`)
    pub beta: [T; 4],
    /// eval₀: folded_eval carried from previous phase
    pub eval_0: [T; 4],
    /// eval₁: sibling evaluation opened from MMCS at this phase
    pub eval_1: [T; 4],

    // Roll-in machinery (extension/base fields)
    /// 0/1 flag indicating whether a reduced opening rolls in at this height
    pub is_roll_in: T,
    /// Reduced opening value for this height (from ReducedOpenings CTL)
    pub reduced_opening: [T; 4],
    /// β² for roll-in computation
    pub beta_squared: [T; 4],
    /// roll_in_contribution = is_roll_in * (β² * reduced_opening)
    pub roll_in_contribution: [T; 4],

    // Index management during folding (base field)
    /// Current leaf index before folding at this phase
    pub domain_index: T,
    /// Parent index after folding (domain_index >> 1)
    pub parent_index: T,
    /// Sibling index at this phase (domain_index ^ 1)
    pub sibling_index: T,
    /// Least-significant bit of domain_index, should be {0,1}
    pub lsb: T,

    // Height (base field)
    /// Current domain height: log_max_height - phase_index - 1
    pub log_height: T,

    // Subgroup points (base field)
    /// x₀ = subgroup_start for this row
    pub x0: T,
    /// x₁ = -x₀
    pub x1: T,

    // FRI folding intermediate arithmetic (extension/base fields)
    /// β - x₀
    pub beta_minus_x0: [T; 4],
    /// eval₁ - eval₀
    pub eval_diff: [T; 4],
    /// x₁ - x₀ = -2 * x₀
    pub x_diff: T,
    /// (x₁ - x₀)⁻¹
    pub x_diff_inv: T,
    /// (β - x₀) * (eval₁ - eval₀)
    pub beta_eval_product: [T; 4],
    /// ((β - x₀)*(eval₁ - eval₀)) * (x₁ - x₀)⁻¹
    pub interpolation_term: [T; 4],
    /// folded_eval = eval₀ + interpolation_term + roll_in_contribution
    pub folded_eval: [T; 4],

    // MMCS commitment for this phase
    /// Commitment (Merkle root digest) to the entire folded layer at this phase
    pub fri_commit: [T; DIGEST_ELEMS],
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
        let (prefix, shorts, suffix) =
            unsafe { self.align_to::<CommitPhaseCols<T, DIGEST_ELEMS>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<T, const DIGEST_ELEMS: usize> BorrowMut<CommitPhaseCols<T, DIGEST_ELEMS>> for [T] {
    fn borrow_mut(&mut self) -> &mut CommitPhaseCols<T, DIGEST_ELEMS> {
        debug_assert_eq!(self.len(), num_commit_phase_cols::<DIGEST_ELEMS>());
        let (prefix, shorts, suffix) =
            unsafe { self.align_to_mut::<CommitPhaseCols<T, DIGEST_ELEMS>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}
