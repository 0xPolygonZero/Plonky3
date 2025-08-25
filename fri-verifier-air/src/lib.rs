//! An AIR implementation for the FRI verifier used in recursive proof systems.

#![no_std]

extern crate alloc;

mod air;
mod columns;
mod generation;

pub use air::*;
pub use columns::*;
pub use generation::*;
