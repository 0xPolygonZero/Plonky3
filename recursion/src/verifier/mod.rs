pub mod circuit_verifier;
pub mod recursive_traits;

use p3_field::BasedVectorSpace;
use p3_uni_stark::{PcsError, StarkGenericConfig, VerificationError, verify as base_verify};

use crate::prover::RecursiveProof;

pub fn verify<SC: StarkGenericConfig, const D: usize>(
    config: &SC,
    proof: RecursiveProof<SC, D>,
) -> Result<(), VerificationError<PcsError<SC>>> {
    assert!(SC::Challenge::DIMENSION == D);
    base_verify(config, &proof.add_air, &proof.add_proof, &vec![])?;
    base_verify(config, &proof.sub_air, &proof.sub_proof, &vec![])
}
