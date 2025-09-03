pub mod tables;

use p3_air::AirBuilder;
use p3_uni_stark::{Proof, StarkGenericConfig, Val};

use crate::air::asic::Asic;
use crate::circuit_builder::gates::event::AllEvents;

pub struct RecursiveProof<SC: StarkGenericConfig> {
    pub proofs: Vec<Proof<SC>>,
}

pub fn prove<SC, AB, const D: usize>(
    config: &SC,
    asic: &Asic<SC, AB, D>,
    all_events: &AllEvents<Val<SC>, D>,
) -> RecursiveProof<SC>
where
    AB: AirBuilder,
    SC: StarkGenericConfig,
{
    let traces = asic.generate_trace(&all_events);

    asic.prove_chips(config, traces)
}

pub trait ProofSystem<const D: usize> {
    type Config: StarkGenericConfig;
    type Builder: AirBuilder<F = Val<Self::Config>>;
    type Proof;
    type Error;

    fn prove<'a>(
        config: &Self::Config,
        asic: &Asic<Self::Config, Self::Builder, D>,
        all_events: &AllEvents<Val<Self::Config>, D>,
    ) -> Result<Self::Proof, Self::Error>;
    fn verify<'a>(
        config: &Self::Config,
        asic: &Asic<Self::Config, Self::Builder, D>,
        proof: &Self::Proof,
    ) -> Result<(), Self::Error>;
}
