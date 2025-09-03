use std::marker::PhantomData;

use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_challenger::DuplexChallenger;
use p3_commit::testing::TrivialPcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_recursion::chips::alu::air::{AddAir, SubAir};
use p3_recursion::chips::asic::Asic;
use p3_recursion::chips::witness::air::WitnessAir;
use p3_recursion::circuit_builder::circuit_builder::{CircuitBuilder, CircuitError};
use p3_recursion::circuit_builder::gates::arith_gates::{AddGate, SubGate};
use p3_recursion::circuit_builder::gates::event::AllEvents;
use p3_recursion::prover::{ProofSystem, RecursiveProof, prove};
use p3_recursion::verifier::verify;
use p3_uni_stark::{DebugConstraintBuilder, PcsError, StarkConfig, Val, VerificationError};
use rand::SeedableRng;
use rand::rngs::SmallRng;

const D: usize = 4;
type Value = BabyBear;
type ExtValue = BinomialExtensionField<Value, 4>;
type Challenge = BinomialExtensionField<Value, 4>;
type Perm = Poseidon2BabyBear<16>;
type Dft = Radix2DitParallel<Value>;
type Challenger = DuplexChallenger<Value, Perm, 16, 8>;
type Pcs = TrivialPcs<Value, Radix2DitParallel<Value>>;
type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;

type SC = StarkConfig<TrivialPcs<Value, Dft>, ExtValue, Challenger>;
type B = DebugConstraintBuilder<'static, Val<SC>>;

pub struct BasicProver {}

impl<const D: usize> ProofSystem<D> for BasicProver {
    type Config = SC;
    type Builder = B;
    type Proof = RecursiveProof<SC>;
    type Error = VerificationError<PcsError<SC>>;

    fn prove<'a>(
        config: &Self::Config,
        asic: &Asic<Self::Config, Self::Builder, D>,
        all_events: &AllEvents<Val<Self::Config>, D>,
    ) -> Result<Self::Proof, Self::Error> {
        let traces = asic.generate_trace(&all_events);

        Ok(asic.prove_chips(config, traces))
    }

    fn verify<'a>(
        config: &Self::Config,
        asic: &Asic<Self::Config, Self::Builder, D>,
        proof: &Self::Proof,
    ) -> Result<(), Self::Error> {
        asic.verify_chips(config, proof)
    }
}

#[test]
pub fn test_simple_circuit() -> Result<(), CircuitError> {
    let mut rng = SmallRng::seed_from_u64(1);
    let perm = Perm::new_from_rng_128(&mut rng);

    let dft = Dft::default();
    let pcs = TrivialPcs {
        dft,
        log_n: 0,
        _phantom: PhantomData,
    };
    let challenger = Challenger::new(perm);
    let config = MyConfig::new(pcs, challenger);

    let mut builder = CircuitBuilder::new();

    let a = builder.new_wire();
    let b = builder.new_wire();
    let c = builder.new_wire();
    let d = builder.new_wire();
    let e = builder.new_wire();

    AddGate::add_to_circuit::<D>(&mut builder, a, b, c);
    SubGate::add_to_circuit::<D>(&mut builder, c, d, e);

    let asic: Asic<MyConfig, DebugConstraintBuilder<Value>, D> = Asic {
        chips: vec![
            Box::new(AddAir::<1>::new()),
            Box::new(SubAir::<1>::new()),
            Box::new(WitnessAir {}),
        ],
    };

    builder.set_wire_value(a, Value::new(10))?;
    builder.set_wire_value(b, Value::new(5))?;
    builder.set_wire_value(d, Value::new(2))?;

    println!("builder wires: {:?}", builder.wires());

    let all_events = builder.generate()?;

    println!("builder wires: {:?}", builder.wires());

    let proof = prove(&config, &asic, &all_events);

    let _proof2 = BasicProver::prove(&config, &asic, &all_events);

    verify(&config, &asic, &proof).unwrap();

    BasicProver::verify(&config, &asic, &proof).unwrap();

    Ok(())
}
