pub mod tables;

use p3_air::AirBuilder;
use p3_uni_stark::{Proof, StarkGenericConfig, Val};

use crate::chips::asic::Asic;
use crate::circuit_builder::gates::event::AllEvents;

pub struct RecursiveProof<SC: StarkGenericConfig> {
    pub proofs: Vec<Proof<SC>>,
}

pub fn prove<SC, AB, const D: usize, const DIGEST_ELEMS: usize>(
    config: &SC,
    asic: &Asic<SC, AB, D, DIGEST_ELEMS>,
    all_events: &AllEvents<Val<SC>, D, DIGEST_ELEMS>,
) -> RecursiveProof<SC>
where
    AB: AirBuilder,
    SC: StarkGenericConfig,
{
    let traces = asic.generate_trace(&all_events);
    println!("trace = ");
    for row in traces[0].row_slices() {
        println!("{:?}", row);
    }

    asic.prove_chips(config, traces)
}

pub trait ProofSystem<const D: usize, const DIGEST_ELEMS: usize> {
    type Config: StarkGenericConfig;
    type Builder: AirBuilder<F = Val<Self::Config>>;
    type Proof;
    type Error;

    fn prove<'a>(
        config: &Self::Config,
        asic: &Asic<Self::Config, Self::Builder, D, DIGEST_ELEMS>,
        all_events: &AllEvents<Val<Self::Config>, D, DIGEST_ELEMS>,
    ) -> Result<Self::Proof, Self::Error>;
    fn verify<'a>(
        config: &Self::Config,
        asic: &Asic<Self::Config, Self::Builder, D, DIGEST_ELEMS>,
        proof: &Self::Proof,
    ) -> Result<(), Self::Error>;
}
