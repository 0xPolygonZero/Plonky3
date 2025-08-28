use std::array;
use std::borrow::Borrow;

use itertools::Itertools;
use p3_air::{Air, AirBuilder, AirBuilderWithPublicValues, BaseAir};
use p3_baby_bear::BabyBear;
use p3_challenger::{HashChallenger, SerializingChallenger32};
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::{BasedVectorSpace, PrimeCharacteristicRing};
use p3_field::{PrimeField64, extension::BinomialExtensionField};
use p3_fri::{HidingFriPcs, create_test_fri_params};
use p3_keccak::{Keccak256Hash, KeccakF};
use p3_matrix::{Matrix, dense::RowMajorMatrix};
use p3_merkle_tree::MerkleTreeHidingMmcs;
use p3_recursion::circuit_builder::{CircuitBuilder, symbolic_to_circuit};
use p3_symmetric::{CompressionFunctionFromHasher, PaddingFreeSponge, SerializingHasher};
use p3_uni_stark::{
    StarkConfig, SymbolicExpression, get_symbolic_constraints, prove, verify,
    verify_with_return_values,
};
use rand::{SeedableRng, rngs::SmallRng};

/// For testing the public values feature
pub struct FibonacciAir {}

const D: usize = 4;
type Val = BabyBear;
type Challenge = BinomialExtensionField<Val, D>;
type Dft = Radix2DitParallel<Val>;

impl<F> BaseAir<F> for FibonacciAir {
    fn width(&self) -> usize {
        NUM_FIBONACCI_COLS
    }
}

impl<AB: AirBuilderWithPublicValues> Air<AB> for FibonacciAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();

        let pis = builder.public_values();

        let a = pis[0];
        let b = pis[1];
        let x = pis[2];

        let (local, next) = (
            main.row_slice(0).expect("Matrix is empty?"),
            main.row_slice(1).expect("Matrix only has 1 row?"),
        );
        let local: &FibonacciRow<AB::Var> = (*local).borrow();
        let next: &FibonacciRow<AB::Var> = (*next).borrow();

        let mut when_first_row = builder.when_first_row();

        when_first_row.assert_eq(local.left.clone(), a);
        when_first_row.assert_eq(local.right.clone(), b);

        let mut when_transition = builder.when_transition();

        // a' <- b
        when_transition.assert_eq(local.right.clone(), next.left.clone());

        // b' <- a + b
        when_transition.assert_eq(local.left.clone() + local.right.clone(), next.right.clone());

        builder.when_last_row().assert_eq(local.right.clone(), x);
    }
}

pub fn generate_trace_rows<F: PrimeField64>(a: u64, b: u64, n: usize) -> RowMajorMatrix<F> {
    assert!(n.is_power_of_two());

    let mut trace = RowMajorMatrix::new(F::zero_vec(n * NUM_FIBONACCI_COLS), NUM_FIBONACCI_COLS);

    let (prefix, rows, suffix) = unsafe { trace.values.align_to_mut::<FibonacciRow<F>>() };
    assert!(prefix.is_empty(), "Alignment should match");
    assert!(suffix.is_empty(), "Alignment should match");
    assert_eq!(rows.len(), n);

    rows[0] = FibonacciRow::new(F::from_u64(a), F::from_u64(b));

    for i in 1..n {
        rows[i].left = rows[i - 1].right;
        rows[i].right = rows[i - 1].left + rows[i - 1].right;
    }

    trace
}

const NUM_FIBONACCI_COLS: usize = 2;

pub struct FibonacciRow<F> {
    pub left: F,
    pub right: F,
}

impl<F> FibonacciRow<F> {
    const fn new(left: F, right: F) -> Self {
        Self { left, right }
    }
}

impl<F> Borrow<FibonacciRow<F>> for [F] {
    fn borrow(&self) -> &FibonacciRow<F> {
        debug_assert_eq!(self.len(), NUM_FIBONACCI_COLS);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<FibonacciRow<F>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

#[test]
fn test_symbolic_to_circuit() {
    type ByteHash = Keccak256Hash;
    let byte_hash = ByteHash {};

    type U64Hash = PaddingFreeSponge<KeccakF, 25, 17, 4>;
    let u64_hash = U64Hash::new(KeccakF {});

    type FieldHash = SerializingHasher<U64Hash>;
    let field_hash = FieldHash::new(u64_hash);

    type MyCompress = CompressionFunctionFromHasher<U64Hash, 2, 4>;
    let compress = MyCompress::new(u64_hash);

    type ValHidingMmcs = MerkleTreeHidingMmcs<
        [Val; p3_keccak::VECTOR_LEN],
        [u64; p3_keccak::VECTOR_LEN],
        FieldHash,
        MyCompress,
        SmallRng,
        4,
        4,
    >;

    let rng = SmallRng::seed_from_u64(1);
    let val_mmcs = ValHidingMmcs::new(field_hash, compress, rng);

    type Challenger = SerializingChallenger32<Val, HashChallenger<u8, ByteHash, 32>>;

    type ChallengeHidingMmcs = ExtensionMmcs<Val, Challenge, ValHidingMmcs>;

    let n = 1 << 3;
    let x = 21;

    let challenge_mmcs = ChallengeHidingMmcs::new(val_mmcs.clone());
    let dft = Dft::default();
    let trace = generate_trace_rows::<Val>(0, 1, n);
    let fri_params = create_test_fri_params(challenge_mmcs, 2);
    type HidingPcs = HidingFriPcs<Val, Dft, ValHidingMmcs, ChallengeHidingMmcs, SmallRng>;
    type MyHidingConfig = StarkConfig<HidingPcs, Challenge, Challenger>;
    let pcs = HidingPcs::new(dft, val_mmcs, fri_params, 4, SmallRng::seed_from_u64(1));
    let challenger = Challenger::from_hasher(vec![], byte_hash);
    let config = MyHidingConfig::new(pcs, challenger);
    let pis = vec![BabyBear::ZERO, BabyBear::ONE, BabyBear::from_u64(x)];
    let proof = prove(&config, &FibonacciAir {}, trace, &pis);
    let _ = verify(&config, &FibonacciAir {}, &proof, &pis).expect("verification failed");

    let air = FibonacciAir {};
    let (alpha, sels, folded_constraints) =
        verify_with_return_values(&config, &air, &proof, &pis).expect("verification failed");

    let symbolic_constraints: Vec<p3_uni_stark::SymbolicExpression<Challenge>> =
        get_symbolic_constraints(&air, 0, pis.len());

    let folded_symbolic_constraints = {
        let mut acc = SymbolicExpression::<Challenge>::Constant(Challenge::ZERO);
        let ch = SymbolicExpression::Constant(alpha);
        for s_c in symbolic_constraints.iter() {
            acc = ch.clone() * acc;
            acc += s_c.clone();
        }
        acc
    };

    let mut circuit = CircuitBuilder::<Val, D>::new();
    let circuit_sels = [
        circuit.new_challenge_wires(),
        circuit.new_challenge_wires(),
        circuit.new_challenge_wires(),
    ];
    let circuit_public_values = [circuit.new_wire(), circuit.new_wire(), circuit.new_wire()];
    let mut circuit_local_values = Vec::with_capacity(NUM_FIBONACCI_COLS);
    let mut circuit_next_values = Vec::with_capacity(NUM_FIBONACCI_COLS);
    for _ in 0..NUM_FIBONACCI_COLS {
        circuit_local_values.push(circuit.new_challenge_wires());
        circuit_next_values.push(circuit.new_challenge_wires());
    }

    // let pis_extension = pis.iter().map(|&x| Challenge::from(x)).collect::<Vec<_>>();
    let sum = symbolic_to_circuit::<Val, Challenge, D>(
        circuit_sels[0].clone(),
        circuit_sels[1].clone(),
        circuit_sels[2].clone(),
        &[],
        &circuit_public_values,
        &[],
        &[],
        &circuit_local_values,
        &circuit_next_values,
        &folded_symbolic_constraints,
        &mut circuit,
    );

    let local_values = &proof.opened_values.trace_local;
    let next_values = &proof.opened_values.trace_next;

    // Set selectors.
    circuit
        .set_challenge_wires(
            circuit_sels[0],
            sels.is_first_row
                .as_basis_coefficients_slice()
                .try_into()
                .unwrap(),
        )
        .unwrap();
    circuit
        .set_challenge_wires(
            circuit_sels[1],
            sels.is_last_row
                .as_basis_coefficients_slice()
                .try_into()
                .unwrap(),
        )
        .unwrap();
    circuit
        .set_challenge_wires(
            circuit_sels[2],
            sels.is_transition
                .as_basis_coefficients_slice()
                .try_into()
                .unwrap(),
        )
        .unwrap();

    // Set public values.
    for (i, pi) in pis.iter().enumerate() {
        circuit
            .set_wire_value(circuit_public_values[i], *pi)
            .unwrap();
    }

    // Set local and next values.
    for (lv, c_lv) in local_values.iter().zip_eq(circuit_local_values) {
        circuit
            .set_challenge_wires(c_lv, lv.as_basis_coefficients_slice().try_into().unwrap())
            .unwrap();
    }
    for (nv, c_nv) in next_values.iter().zip_eq(circuit_next_values) {
        circuit
            .set_challenge_wires(c_nv, nv.as_basis_coefficients_slice().try_into().unwrap())
            .unwrap();
    }

    let _ = circuit.generate();

    let challenge_sum = Challenge::from_basis_coefficients_slice(&array::from_fn::<_, D, _>(|i| {
        circuit.get_wire_value(sum[i]).unwrap().unwrap()
    }))
    .unwrap();

    assert!(folded_constraints == challenge_sum);
}
