use std::ops::{Add, Div, Mul, Neg};

/// Mimic some features of a ring, without ring axioms (Abelian group, Associativity, Distributivity)
pub trait Numeric:
    Sized + AddIdentity + Add<Output = Self> + Mul<Output = Self> + Neg<Output = Self>
{
}
// TODO: rename to make sure numeric is not meant to be math numeric

pub trait AddIdentity: Sized {
    fn zero() -> Self;
}

pub trait NumericAddIdentity: AddIdentity {
    /// For numerical reasons
    fn is_zero(&self, numerical_threshold: Option<f64>) -> bool;
}

pub trait MulIdentity: Sized {
    fn one() -> Self;
}

pub trait NumericMulIdentity: MulIdentity {
    /// For numerical reasons
    fn is_one(&self, numerical_threshold: Option<f64>) -> bool;
}

pub trait NumericInvertible {
    fn is_invertible(&self) -> bool;
}

pub trait Inverse: Sized {
    type Error;
    fn inverse(&self) -> std::result::Result<Self, Self::Error>;
}

pub trait SquareRoot {
    fn square_root(&self) -> Self;
}

pub trait Signed {
    /// Cannot call it `signum` due to naming conflict.
    fn signed(&self) -> Self;
}

pub trait IsPositive {
    /// Cannot call it `signum` due to naming conflict.
    fn is_positive(&self) -> bool;
}

pub trait MuliplicativeNumeric:
    Sized + MulIdentity + Mul<Output = Self> + Div<Output = Self>
{
} // Inverse {}

#[macro_export]
macro_rules! impl_add_identity {
    ($impl_type:ty) => {
        impl AddIdentity for $impl_type {
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
        impl MulIdentity for $impl_type {
            fn one() -> Self {
                1 as $impl_type
            }
        }
    };
}

#[macro_export]
macro_rules! impl_square_root {
    ($impl_type:ty) => {
        impl SquareRoot for $impl_type {
            fn square_root(&self) -> Self {
                self.abs().sqrt()
            }
        }
    };
}

#[macro_export]
macro_rules! impl_signed {
    ($impl_type:ty) => {
        impl Signed for $impl_type {
            fn signed(&self) -> Self {
                self.signum()
            }
        }
    };
}

#[macro_export]
macro_rules! impl_is_positive {
    ($impl_type:ty) => {
        impl IsPositive for $impl_type {
            fn is_positive(&self) -> bool {
                self.signum() > (0 as $impl_type)
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

// implement MuliplicativeNumeric
impl MuliplicativeNumeric for f32 {}
impl MuliplicativeNumeric for f64 {}

// implement SignedSquareRoot
impl_square_root! { f32 }
impl_square_root! { f64 }

// implement Signed
impl_signed! { f32 }
impl_signed! { f64 }

// implement IsPositive
impl_is_positive! { i8 }
impl_is_positive! { i16 }
impl_is_positive! { i32 }
impl_is_positive! { i64 }
impl_is_positive! { f32 }
impl_is_positive! { f64 }

// TODO: add bigdecimal support
// TODO: add 'complex number' support
