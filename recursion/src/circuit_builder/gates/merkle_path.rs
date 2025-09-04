use std::char::MAX;

use itertools::Itertools;
use p3_field::{Field, PrimeField, PrimeField32};
use p3_symmetric::FieldCompression;

use crate::chips::merkle_path::air::{MerklePathEvent, SiblingWithExtraOpening, number_of_nodes};
use crate::circuit_builder::gates::event::AllEvents;
use crate::circuit_builder::gates::gate::Gate;
use crate::circuit_builder::{CircuitBuilder, CircuitError, WireId};

pub struct MerklePathGate<F, C, const DIGEST_ELEMS: usize, const MAX_TREE_HEIGHT: usize>
where
    F: PrimeField32,
    C: FieldCompression<F, 2, DIGEST_ELEMS> + 'static,
{
    leaf: [WireId; DIGEST_ELEMS],
    index: WireId,
    siblings: Vec<SiblingWithExtraOpening<[WireId; DIGEST_ELEMS]>>,
    root: [WireId; DIGEST_ELEMS],
    compress: C,
    _marker: std::marker::PhantomData<F>,
}

impl<F, C, const DIGEST_ELEMS: usize, const MAX_TREE_HEIGHT: usize>
    MerklePathGate<F, C, DIGEST_ELEMS, MAX_TREE_HEIGHT>
where
    F: PrimeField32,
    C: FieldCompression<F, 2, DIGEST_ELEMS>,
{
    pub fn new(
        leaf: [WireId; DIGEST_ELEMS],
        index: WireId,
        siblings: Vec<SiblingWithExtraOpening<[WireId; DIGEST_ELEMS]>>,
        root: [WireId; DIGEST_ELEMS],
        compress: C,
    ) -> Self {
        println!("sbl len = {:?}, MTH  = {MAX_TREE_HEIGHT}", siblings.len());
        debug_assert!(siblings.len() <= MAX_TREE_HEIGHT);

        MerklePathGate {
            leaf,
            index,
            siblings,
            root,
            compress,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn add_to_circuit<const D: usize>(
        builder: &mut CircuitBuilder<F, D, DIGEST_ELEMS>,
        leaf: [WireId; DIGEST_ELEMS],
        index: WireId,
        siblings: Vec<SiblingWithExtraOpening<[WireId; DIGEST_ELEMS]>>,
        root: [WireId; DIGEST_ELEMS],
        compress: C,
    ) {
        let gate = Self::new(leaf, index, siblings, root, compress);
        builder.add_gate(Box::new(gate));
    }

    pub fn get_root(
        compress: &C,
        leaf: &[F; DIGEST_ELEMS],
        index: &u32,
        siblings: &Vec<SiblingWithExtraOpening<[F; DIGEST_ELEMS]>>,
    ) -> [F; DIGEST_ELEMS] {
        siblings
            .iter()
            .zip((0..32).map(|i| (index >> i & 1) == 1))
            .fold(*leaf, |state, ((sibling, extra_row), index_bit)| {
                let mut state = if index_bit {
                    compress.compress_field([*sibling, state])
                } else {
                    compress.compress_field([state, *sibling])
                };
                if let Some(extra_row) = extra_row {
                    state = compress.compress_field([state, *extra_row])
                }
                state
            })
    }
}

impl<F, C, const D: usize, const DIGEST_ELEMS: usize, const MAX_TREE_HEIGHT: usize>
    Gate<F, D, DIGEST_ELEMS> for MerklePathGate<F, C, DIGEST_ELEMS, MAX_TREE_HEIGHT>
where
    F: PrimeField32,
    C: FieldCompression<F, 2, DIGEST_ELEMS>,
{
    fn n_inputs(&self) -> usize {
        self.leaf.len() + DIGEST_ELEMS * self.siblings.iter().map(number_of_nodes).sum::<usize>()
    }

    fn n_outputs(&self) -> usize {
        self.root.len()
    }

    fn generate(
        &self,
        builder: &mut CircuitBuilder<F, D, DIGEST_ELEMS>,
        all_events: &mut AllEvents<F, D, DIGEST_ELEMS>,
    ) -> Result<(), CircuitError> {
        println!("gen merkle path");
        let get_wire_value = |wire_name| {
            let builder_ref = &builder;
            let msg = format!("Wire not set: MerklePathGate::{:?}", wire_name);
            move |wire| {
                builder_ref.get_wire_value(wire).expect(&msg).expect(&msg) // TOD: Manage errors  properly
            }
        };
        let leaf = self.leaf.map(get_wire_value("leaf"));
        let index = get_wire_value("index")(self.index);
        let siblings = self
            .siblings
            .iter()
            .map(|(sibling_wires, extra_row_wires)| {
                (
                    sibling_wires.map(get_wire_value("sibling")),
                    extra_row_wires
                        .map(|extra_row_wires| extra_row_wires.map(get_wire_value("extra row"))),
                )
            })
            .collect_vec();
        let index = index.as_canonical_u32();
        let root = Self::get_root(&self.compress, &leaf, &index, &siblings);
        builder.set_wire_values(&self.root, &root)?;

        println!(
            "get_root(compress, {:?}, {:?}, {:?}) = {:?}",
            leaf, index, siblings, root
        );

        all_events.merkle_path_events.push(MerklePathEvent {
            leaf,
            index,
            siblings,
            root,
        });
        println!("merkle path generated");
        Ok(())
    }
}

#[cfg(test)]
mod test {

    use p3_baby_bear::BabyBear;
    use p3_challenger::{HashChallenger, SerializingChallenger32};
    use p3_commit::ExtensionMmcs;
    use p3_dft::Radix2DitParallel;
    use p3_field::PrimeCharacteristicRing;
    use p3_field::extension::BinomialExtensionField;
    use p3_fri::{TwoAdicFriPcs, create_benchmark_fri_params};
    use p3_keccak::{Keccak256Hash, KeccakF};
    use p3_merkle_tree::MerkleTreeMmcs;
    use p3_symmetric::{
        CompressionFunctionFromHasher, FieldCompression, PaddingFreeSponge, SerializingHasher,
    };
    use p3_uni_stark::{DebugConstraintBuilder, StarkConfig};
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    use crate::chips::alu::air::AddAir;
    use crate::chips::asic::Asic;
    use crate::chips::merkle_path::air::{
        MerklePathAir, SiblingWithExtraOpening, make_counter, make_field_counter,
    };
    use crate::chips::witness::air::WitnessAir;
    use crate::circuit_builder::gates::arith_gates::AddGate;
    use crate::circuit_builder::gates::merkle_path::MerklePathGate;
    use crate::circuit_builder::{CircuitBuilder, CircuitError, WireId};
    use crate::prover::prove;
    use crate::verifier::verify;

    const D: usize = 4;
    const DIGEST_ELEMS: usize = 4;
    const MAX_TREE_HEIGHT: usize = 8;
    type Value = BabyBear;
    type Challenge = BinomialExtensionField<Value, 4>;
    type ByteHash = Keccak256Hash;
    type Dft = Radix2DitParallel<Value>;
    type MmcsHash = PaddingFreeSponge<KeccakF, 25, 17, 4>;
    type FieldHash = SerializingHasher<MmcsHash>;
    type MmcsCompress = CompressionFunctionFromHasher<MmcsHash, 2, 4>;
    type ValMmcs = MerkleTreeMmcs<
        [Value; p3_keccak::VECTOR_LEN],
        [u64; p3_keccak::VECTOR_LEN],
        FieldHash,
        MmcsCompress,
        4,
    >;
    type ChallengeMmcs = ExtensionMmcs<Value, Challenge, ValMmcs>;
    type Challenger = SerializingChallenger32<Value, HashChallenger<u8, ByteHash, 32>>;
    type Pcs = TwoAdicFriPcs<Value, Dft, ValMmcs, ChallengeMmcs>;
    type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;
    type U64Hash = PaddingFreeSponge<KeccakF, 25, 17, DIGEST_ELEMS>;
    type MyCompress = CompressionFunctionFromHasher<U64Hash, 2, DIGEST_ELEMS>;

    #[test]
    pub fn test_merkle_path_gate() -> Result<(), CircuitError> {
        let (config, compress) = get_config_and_compress();

        let mut builder = CircuitBuilder::new();

        let siblings: Vec<SiblingWithExtraOpening<[WireId; DIGEST_ELEMS]>> = (0..MAX_TREE_HEIGHT)
            .map(|i| {
                (
                    std::array::from_fn(|_| builder.new_wire()),
                    if i == MAX_TREE_HEIGHT / 2 - 1 {
                        Some(std::array::from_fn(|_| builder.new_wire()))
                    } else {
                        None
                    },
                )
            })
            .collect();
        let leaf: [WireId; DIGEST_ELEMS] = std::array::from_fn(|_| builder.new_wire());
        let index = builder.new_wire();
        let root: [WireId; DIGEST_ELEMS] = std::array::from_fn(|_| builder.new_wire());

        MerklePathGate::<Value, _, DIGEST_ELEMS, MAX_TREE_HEIGHT>::add_to_circuit(
            &mut builder,
            leaf,
            index,
            siblings.clone(),
            root,
            compress.clone(),
        );

        let asic: Asic<MyConfig, DebugConstraintBuilder<Value>, D, _> = Asic {
            chips: vec![
                Box::new(MerklePathAir::<_, _, _, MAX_TREE_HEIGHT>::new(compress)),
                Box::new(WitnessAir {}),
            ],
        };

        let mut field_counter = make_field_counter::<_, 4>();
        let mut index_counter = make_counter();

        // for &wire in leaf.iter().chain(std::iter::once(&index)) {
        //     builder.set_wire_value(wire, rng.random())?;
        // }
        builder.set_wire_values(&leaf, &field_counter())?;
        builder.set_wire_value(index, Value::from_u32(index_counter()))?;

        for (wires, extra_wires) in siblings.iter() {
            // builder.set_wire_values(wires, &std::array::from_fn::<_, 4, _>(|_| rng.random()))?;
            builder.set_wire_values(wires, &field_counter())?;
            extra_wires.map(|extra_wires| {
                builder.set_wire_values(
                    &extra_wires,
                    // &std::array::from_fn::<_, 4, _>(|_| rng.random()),
                    &field_counter(),
                )
            });
        }

        let all_events = builder.generate()?;

        println!("all_events = {:?}", all_events.merkle_path_events);

        let proof = prove(&config, &asic, &all_events);

        verify(&config, &asic, &proof).unwrap();

        Ok(())
    }

    fn get_config_and_compress() -> (MyConfig, MyCompress) {
        let mut _rng = SmallRng::seed_from_u64(1);
        let byte_hash = ByteHash {};

        let u64_hash = U64Hash::new(KeccakF {});

        let field_hash = FieldHash::new(u64_hash);
        let compress = MyCompress::new(u64_hash);
        let val_mmcs = ValMmcs::new(field_hash, compress);
        let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());
        let fri_params = create_benchmark_fri_params(challenge_mmcs);

        let dft = Dft::default();
        let pcs = Pcs::new(dft, val_mmcs, fri_params);
        let challenger = Challenger::from_hasher(vec![], byte_hash);
        let config = MyConfig::new(pcs, challenger);

        let u64_hash = U64Hash::new(KeccakF {});

        let compress = MyCompress::new(u64_hash);

        println!(
            "compres_filed(0, 0) = {:?}",
            compress.compress_field([[Value::ZERO; 4], [Value::ZERO; 4]])
        );
        (config, compress)
    }

    pub const ANOTHER_MAX_TREE_HEIGHT: usize = 24;
    #[test]
    pub fn test_two_merkle_paths() -> Result<(), CircuitError> {
        let (config, compress) = get_config_and_compress();

        let mut builder = CircuitBuilder::new();

        let siblings: Vec<SiblingWithExtraOpening<[WireId; DIGEST_ELEMS]>> = (0..MAX_TREE_HEIGHT)
            .map(|i| {
                (
                    std::array::from_fn(|_| builder.new_wire()),
                    if i == MAX_TREE_HEIGHT / 2 - 1 {
                        Some(std::array::from_fn(|_| builder.new_wire()))
                    } else {
                        None
                    },
                )
            })
            .collect();
        let leaf: [WireId; DIGEST_ELEMS] = std::array::from_fn(|_| builder.new_wire());
        let index = builder.new_wire();
        let root: [WireId; DIGEST_ELEMS] = std::array::from_fn(|_| builder.new_wire());

        // Add a new gate checking a merkle path of size 8.
        MerklePathGate::<_, _, DIGEST_ELEMS, MAX_TREE_HEIGHT>::add_to_circuit::<D>(
            &mut builder,
            leaf,
            index,
            siblings.clone(),
            root.clone(),
            compress.clone(),
        );

        let another_siblings: Vec<SiblingWithExtraOpening<[WireId; DIGEST_ELEMS]>> = (0
            ..ANOTHER_MAX_TREE_HEIGHT)
            .map(|i| {
                (
                    std::array::from_fn(|_| builder.new_wire()),
                    if i == ANOTHER_MAX_TREE_HEIGHT / 2 - 1 {
                        Some(std::array::from_fn(|_| builder.new_wire()))
                    } else {
                        None
                    },
                )
            })
            .collect();
        let another_leaf: [WireId; DIGEST_ELEMS] = root;
        let another_index = builder.new_wire();
        let another_root: [WireId; DIGEST_ELEMS] = std::array::from_fn(|_| builder.new_wire());

        // The other index is add_index + 1
        let one = builder.add_constant(Value::ONE);
        AddGate::add_to_circuit(&mut builder, index, one, another_index);

        // Add another gate cheking a merkle path of size 24 with leaf the root of the previous gate
        MerklePathGate::<Value, _, DIGEST_ELEMS, ANOTHER_MAX_TREE_HEIGHT>::add_to_circuit(
            &mut builder,
            another_leaf,
            another_index,
            another_siblings.clone(),
            another_root,
            compress.clone(),
        );

        let mut field_counter = make_field_counter::<_, 4>();
        let mut index_counter = make_counter();

        builder.set_wire_values(&leaf, &field_counter())?;
        builder.set_wire_value(index, Value::from_u32(index_counter()))?;

        for (wires, extra_wires) in siblings.iter() {
            builder.set_wire_values(wires, &field_counter())?;
            extra_wires.map(|extra_wires| builder.set_wire_values(&extra_wires, &field_counter()));
        }

        for (wires, extra_wires) in another_siblings.iter() {
            builder.set_wire_values(wires, &field_counter())?;
            extra_wires.map(|extra_wires| builder.set_wire_values(&extra_wires, &field_counter()));
        }

        let asic: Asic<MyConfig, DebugConstraintBuilder<Value>, D, _> = Asic {
            chips: vec![
                Box::new(MerklePathAir::<_, _, _, ANOTHER_MAX_TREE_HEIGHT>::new(
                    compress,
                )),
                Box::new(WitnessAir {}),
                Box::new(AddAir::<1>::new()),
            ],
        };

        let all_events = builder.generate()?;

        println!("all_events = {:?}", all_events.merkle_path_events);

        let proof = prove(&config, &asic, &all_events);

        verify(&config, &asic, &proof).unwrap();

        Ok(())
    }
}
