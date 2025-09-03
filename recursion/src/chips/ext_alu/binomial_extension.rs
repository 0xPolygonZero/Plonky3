use std::borrow::{Borrow, BorrowMut};
use std::ops::{Add, Mul, Sub};

use p3_field::{Algebra, Field};
use rand::Rng;
use rand::distr::{Distribution, StandardUniform};

/// A binomial extension element represented over a generic type `T`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(C)]
pub struct BinomialExtension<T, const D: usize>(pub [T; D]);

impl<F, const D: usize> BinomialExtension<F, D> {
    /// Creates a new binomial extension element from a base element.
    pub fn from_base(b: F) -> Self
    where
        F: Field,
    {
        let mut arr: [F; D] = core::array::from_fn(|_| F::ZERO);
        arr[0] = b;
        Self(arr)
    }

    /// Returns a reference to the underlying slice.
    pub const fn as_base_slice(&self) -> &[F] {
        &self.0
    }

    /// Creates a new binomial extension element from a binomial extension element.
    pub fn from<S: Into<F> + Clone>(from: &BinomialExtension<S, D>) -> Self {
        BinomialExtension(core::array::from_fn(|i| from.0[i].clone().into()))
    }
}

impl<F, const D: usize> Borrow<BinomialExtension<F, D>> for [F] {
    fn borrow(&self) -> &BinomialExtension<F, D> {
        debug_assert_eq!(self.len(), D);
        let (prefix, shorts, _suffix) = unsafe { self.align_to::<BinomialExtension<F, D>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<F, const D: usize> BorrowMut<BinomialExtension<F, D>> for [F] {
    fn borrow_mut(&mut self) -> &mut BinomialExtension<F, D> {
        debug_assert_eq!(self.len(), D);
        let (prefix, shorts, _suffix) = unsafe { self.align_to_mut::<BinomialExtension<F, D>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &mut shorts[0]
    }
}

impl<T: Add<Output = T> + Clone, const D: usize> Add for BinomialExtension<T, D> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(core::array::from_fn(|i| {
            self.0[i].clone() + rhs.0[i].clone()
        }))
    }
}

impl<T: Sub<Output = T> + Clone, const D: usize> Sub for BinomialExtension<T, D> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(core::array::from_fn(|i| {
            self.0[i].clone() - rhs.0[i].clone()
        }))
    }
}

impl<F: Add<Output = F> + Mul<Output = F> + Algebra<F>, const D: usize> Mul
    for BinomialExtension<F, D>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = [F::ZERO; D];
        // This value is specific for BabyBear prime's extension `F_p[x]/(x^4 - 11)`.
        let w = F::from_u8(11); // TODO: Is this correct?

        for i in 0..D {
            for j in 0..D {
                if i + j >= D {
                    result[i + j - D] = result[i + j - D].clone()
                        + w.clone() * self.0[i].clone() * rhs.0[j].clone();
                } else {
                    result[i + j] = result[i + j].clone() + self.0[i].clone() * rhs.0[j].clone();
                }
            }
        }

        Self(result)
    }
}

impl<F, const D: usize> Distribution<BinomialExtension<F, D>> for StandardUniform
where
    StandardUniform: Distribution<F>,
{
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> BinomialExtension<F, D> {
        BinomialExtension(core::array::from_fn(|_| rng.random()))
    }
}
