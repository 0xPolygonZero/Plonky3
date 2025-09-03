use p3_air::AirBuilder;
use p3_uni_stark::{PcsError, StarkGenericConfig, VerificationError};

use crate::air::asic::Asic;
use crate::prover::RecursiveProof;

pub mod circuit_verifier;
pub mod recursive_traits;

pub fn verify<AB: AirBuilder, SC: StarkGenericConfig, const D: usize>(
    config: &SC,
    asic: &Asic<SC, AB, D>,
    proof: &RecursiveProof<SC>,
) -> Result<(), VerificationError<PcsError<SC>>> {
    asic.verify_chips(config, proof)
}
