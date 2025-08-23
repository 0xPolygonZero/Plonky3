//! An AIR for the FRI verifier.

#![no_std]

extern crate alloc;

mod air;
mod columns;
mod generation;

pub use air::*;
pub use columns::*;
pub use generation::*;

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use p3_air::BaseAir;
    use p3_field::PrimeCharacteristicRing;
    use p3_matrix::Matrix;

    type TestField = p3_baby_bear::BabyBear;

    #[test]
    fn test_commit_phase_step_creation() {
        let step = create_commit_phase_step::<TestField>(
            0,                        // query_id
            0,                        // phase_idx
            TestField::from_u32(42),  // beta
            8,                        // domain_index
            TestField::from_u32(100), // eval_0
            TestField::from_u32(200), // eval_1
            TestField::from_u32(7),   // x0
            TestField::ZERO,          // roll_in_value
            4,                        // log_height
            true,                     // mmcs_verified
        );

        // Verify computed fields
        assert_eq!(step.parent_index, TestField::from_usize(4)); // 8 >> 1
        assert_eq!(step.sibling_index, TestField::from_usize(9)); // 8 ^ 1
        assert_eq!(step.x1, -TestField::from_u32(7)); // -x0
    }

    #[test]
    fn test_trace_generation() {
        let steps = vec![
            create_commit_phase_step::<TestField>(
                0,
                0,
                TestField::from_u32(42),
                8,
                TestField::from_u32(100),
                TestField::from_u32(200),
                TestField::from_u32(7),
                TestField::ZERO,
                4,
                true,
            ),
            create_commit_phase_step::<TestField>(
                0,
                1,
                TestField::from_u32(43),
                4,
                TestField::from_u32(150),
                TestField::from_u32(250),
                TestField::from_u32(3),
                TestField::ZERO,
                3,
                true,
            ),
        ];

        let trace = generate_commit_phase_trace(steps);
        assert_eq!(trace.height(), 2);
        assert_eq!(trace.width(), core::mem::size_of::<CommitPhaseCols<u8>>());
    }

    #[test]
    fn test_air_width_consistency() {
        let air = CommitPhaseAir::<TestField>::new();
        let expected_width = core::mem::size_of::<CommitPhaseCols<u8>>();

        assert_eq!(air.width(), expected_width);
        assert_eq!(air.width(), num_commit_phase_cols());
    }

    #[test]
    fn test_commit_phase_cols_size() {
        // Verify the column count matches the struct size
        assert_eq!(num_commit_phase_cols(), 20); // Should be 20 fields
        assert_eq!(core::mem::size_of::<CommitPhaseCols<u8>>(), 20);
    }
}
