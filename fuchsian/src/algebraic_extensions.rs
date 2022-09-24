use std::ops::{Add, Mul, Neg};

/// Mimic some features of a ring, without ring axioms (Abelian group, Associativity, Distributivity)
pub trait Numeric:
    Sized + AddIdentityElement + Add<Output = Self> + Mul<Output = Self> + Neg<Output = Self>
{
}
// TODO: rename to make sure numeric is not meant to be math numeric

pub trait AddIdentityElement: Sized {
    fn zero() -> Self;
}

pub trait NumericAddIdentity: AddIdentityElement {
    /// For numerical reasons
    fn is_zero(&self, threshold: Option<f64>) -> bool;
}

pub trait MulIdentityElement: Sized {
    fn one() -> Self;
}

pub trait NumericMulIdentity: MulIdentityElement {
    /// For numerical reasons
    fn is_one(&self, threshold: Option<f64>) -> bool;
}

pub trait NumericInvertible {
    fn is_invertible(&self) -> bool;
}

// impl<T> NumericMulIdentity for T
// where
//     T: MulIdentity + NumericAddIdentity + Add<Output = Self> + Neg<Output = Self> + Copy,
// {
//     fn is_one(&self, threshold: Option<f64>) -> bool {
//         (*self + (-T::one())).is_zero(threshold)
//     }
// }

#[macro_export]
macro_rules! impl_add_identity {
    ($impl_type:ty) => {
        impl AddIdentityElement for $impl_type {
            fn zero() -> Self {
                0 as $impl_type
            }
        }
    };
}

#[macro_export]
macro_rules! impl_numeric_add_identity {
    ($impl_type:ty) => {
        impl NumericAddIdentity for $impl_type {
            fn is_zero(&self, threshold: Option<f64>) -> bool {
                match threshold {
                    Some(tol) => self.abs() as f64 <= tol,
                    None => self.abs() == Self::zero(),
                }
            }
        }
    };
}

#[macro_export]
macro_rules! impl_numeric_mul_identity {
    ($impl_type:ty) => {
        impl NumericMulIdentity for $impl_type {
            fn is_one(&self, threshold: Option<f64>) -> bool {
                (*self - Self::one()).is_zero(threshold)
            }
        }
    };
}

#[macro_export]
macro_rules! impl_mul_identity {
    ($impl_type:ty) => {
        impl MulIdentityElement for $impl_type {
            fn one() -> Self {
                1 as $impl_type
            }
        }
    };
}

// implement AddIdentity
impl_add_identity! { i8 }
impl_add_identity! { i16 }
impl_add_identity! { i32 }
impl_add_identity! { i64 }
impl_add_identity! { f32 }
impl_add_identity! { f64 }

// implement NumericAddIdentity
impl_numeric_add_identity! { i8 }
impl_numeric_add_identity! { i16 }
impl_numeric_add_identity! { i32 }
impl_numeric_add_identity! { i64 }
impl_numeric_add_identity! { f32 }
impl_numeric_add_identity! { f64 }

// implement MulIdentity
impl_mul_identity! { i8 }
impl_mul_identity! { i16 }
impl_mul_identity! { i32 }
impl_mul_identity! { i64 }
impl_mul_identity! { f32 }
impl_mul_identity! { f64 }

// implement NumericAddIdentity
impl_numeric_mul_identity! { i8 }
impl_numeric_mul_identity! { i16 }
impl_numeric_mul_identity! { i32 }
impl_numeric_mul_identity! { i64 }
impl_numeric_mul_identity! { f32 }
impl_numeric_mul_identity! { f64 }

// implement Numeric
impl Numeric for i8 {}
impl Numeric for i16 {}
impl Numeric for i32 {}
impl Numeric for i64 {}
impl Numeric for f32 {}
impl Numeric for f64 {}

// TODO: add bigdecimal support
// TODO: add 'complex number' support
