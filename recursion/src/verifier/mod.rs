use p3_air::AirBuilder;
use p3_uni_stark::{PcsError, StarkGenericConfig, VerificationError};

use crate::chips::asic::Asic;
use crate::prover::RecursiveProof;

pub mod circuit_verifier;
pub mod recursive_pcs;
pub mod recursive_traits;

pub fn verify<AB: AirBuilder, SC: StarkGenericConfig, const D: usize, const DIGEST_ELEMS: usize>(
    config: &SC,
    asic: &Asic<SC, AB, D, DIGEST_ELEMS>,
    proof: &RecursiveProof<SC>,
) -> Result<(), VerificationError<PcsError<SC>>> {
    asic.verify_chips(config, proof)
}
