use crate::algebraic_extensions::{
    AddIdentity, Inverse, MulIdentity, Numeric, NumericAddIdentity, NumericMulIdentity,
};
use std::{
    fmt,
    ops::{Add, Div, Mul, Neg, Sub},
};

// TODO: consider https://docs.rs/katex-doc/latest/katex_doc/ or
//  https://crates.io/crates/rust-latex-doc-minimal-example for math formulas

/// https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
/// Corresponds to the 2x2 matrix-vector multiplication of the matrix
/// $[a, b; c, d] * [x; y]$
/// for $[x; y]$ in the (real) Euclidean vector space.
/// Corresponds to the complex-valued function
/// $z complex -> f(z) = \frac{a*z + b}{cz + d}$
pub struct MoebiusTransformation<T> {
    pub a: T,
    pub b: T,
    pub c: T,
    pub d: T,
}

#[macro_export]
macro_rules! moebius {
    ($impl_type:ty, $a:expr, $b:expr, $c:expr, $d:expr) => {
        MoebiusTransformation::<$impl_type>::new($a, $b, $c, $d)
    };
}

impl<T> MoebiusTransformation<T> {
    pub fn new(a: T, b: T, c: T, d: T) -> Self {
        Self { a, b, c, d }
    }

    /// The square root of the `squared_sum` is its length within the vector space of 2x2 matrices.
    pub fn squared_sum(&self) -> T
    where
        T: Numeric + std::marker::Copy,
    {
        self.a * self.a + self.b * self.b + self.c * self.c + self.d * self.d
    }

    pub fn determinant(&self) -> T
    where
        T: Numeric + std::marker::Copy,
    {
        self.a * self.d + (-(self.b * self.c))
    }

    pub fn is_invertible(&self, threshold: Option<f64>) -> bool
    where
        T: Numeric + NumericAddIdentity + std::marker::Copy,
    {
        !self.determinant().is_zero(threshold)
    }

    pub fn inverse(&self, threshold: Option<f64>) -> Option<Self>
    where
        T: Numeric + NumericAddIdentity + std::marker::Copy + Div<Output = T>,
    {
        if self.is_invertible(threshold) {
            // TODO: could calc the determinant only once (already done in is_invertible)
            let det = self.determinant();
            return Some(Self {
                a: self.d / det,
                b: -self.b / det,
                c: -self.c / det,
                d: self.a / det,
            });
        }
        None
    }
}

// ########################
// Implement useful traits
// ########################

impl<T> fmt::Display for MoebiusTransformation<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MT[{}, {}; {}, {}]", self.a, self.b, self.c, self.d)
    }
}

impl<T> Copy for MoebiusTransformation<T> where T: Copy {}
impl<T> Clone for MoebiusTransformation<T>
where
    T: Clone,
{
    fn clone(&self) -> Self {
        Self::new(
            self.a.clone(),
            self.b.clone(),
            self.c.clone(),
            self.d.clone(),
        )
    }
}

impl<T> PartialEq for MoebiusTransformation<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b && self.c == other.c && self.d == other.d
    }
}
impl<T> Eq for MoebiusTransformation<T> where T: PartialEq {}

impl<T> std::hash::Hash for MoebiusTransformation<T>
where
    T: std::hash::Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.a.hash(state);
        self.b.hash(state);
        self.c.hash(state);
        self.d.hash(state);
    }
}

// ########################
// Establish algebraic structures on the set of Moebius-transformations (group, ring, vector space)
// ########################

impl<T> AddIdentity for MoebiusTransformation<T>
where
    T: AddIdentity,
{
    fn zero() -> Self {
        MoebiusTransformation::new(T::zero(), T::zero(), T::zero(), T::zero())
    }
}

impl<T> NumericAddIdentity for MoebiusTransformation<T>
where
    T: NumericAddIdentity,
{
    fn is_zero(&self, threshold: Option<f64>) -> bool {
        self.a.is_zero(threshold)
            && self.b.is_zero(threshold)
            && self.c.is_zero(threshold)
            && self.d.is_zero(threshold)
    }
}

/// Corresponds to the matrix [1, 0; 0, 1]
impl<T> MulIdentity for MoebiusTransformation<T>
where
    T: MulIdentity + AddIdentity,
{
    fn one() -> Self {
        MoebiusTransformation::new(T::one(), T::zero(), T::zero(), T::one())
    }
}

impl<T> NumericMulIdentity for MoebiusTransformation<T>
where
    T: NumericMulIdentity + NumericAddIdentity,
{
    fn is_one(&self, threshold: Option<f64>) -> bool {
        self.a.is_one(threshold)
            && self.b.is_zero(threshold)
            && self.c.is_zero(threshold)
            && self.d.is_one(threshold)
    }
}

impl<T> Add for MoebiusTransformation<T>
where
    T: Numeric,
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            a: self.a + other.a,
            b: self.b + other.b,
            c: self.c + other.c,
            d: self.d + other.d,
        }
    }
}

impl<T> Neg for MoebiusTransformation<T>
where
    T: Numeric + Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            a: -self.a,
            b: -self.b,
            c: -self.c,
            d: -self.d,
        }
    }
}

impl<T> Sub for MoebiusTransformation<T>
where
    T: Numeric + Neg<Output = T>,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl<T> Mul for MoebiusTransformation<T>
where
    T: Numeric + Copy,
{
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            a: self.a * other.a + self.b * other.c,
            b: self.a * other.b + self.b * other.d,
            c: self.c * other.a + self.d * other.c,
            d: self.c * other.b + self.d * other.d,
        }
    }
}

// Scalability, corresponds to scalar * matrix
impl<T> Mul<T> for MoebiusTransformation<T>
where
    T: Mul<Output = T> + Copy,
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Self {
            a: rhs * self.a,
            b: rhs * self.b,
            c: rhs * self.c,
            d: rhs * self.d,
        }
    }
}

// TODO: move
pub const DEFAULT_THRESHOLD: f64 = 1e-16;

impl<T> Inverse for MoebiusTransformation<T>
where
    T: Numeric + NumericAddIdentity + std::marker::Copy + Div<Output = T>,
{
    type Error = &'static str;

    fn inverse(&self) -> std::result::Result<Self, Self::Error> {
        if let Some(m) = self.inverse(Some(DEFAULT_THRESHOLD)) {
            return Ok(m);
        }
        Err("Moebius transformation is not invertible. Determinant smaller than 1e-16")
    }
}

#[cfg(test)]
mod tests {
    use super::MoebiusTransformation;
    use crate::algebraic_extensions::{AddIdentity, MulIdentity, NumericMulIdentity};

    #[test]
    fn test_macro() {
        let m = moebius!(i8, 1, 2, 3, 4);
        assert_eq!(m.a, 1);
        assert_eq!(m.d, 4);

        let m = moebius!(i64, 1, 2, 3, 4);
        assert_eq!(m.a, 1);
        assert_eq!(m.d, 4);

        let m = moebius!(f32, 1.0, 2.0, 3.0, 4.0);
        assert_eq!(m.a, 1.0);
        assert_eq!(m.d, 4.0);
    }

    #[test]
    fn test_is_invertible() {
        let m1 = moebius!(f32, 1.0, 2.0, 3.0, 4.0);
        assert_eq!(m1.determinant(), -2.0);
        assert!(m1.is_invertible(None));

        let m2 = moebius!(f64, 1.0, 2.0, 2.0, 4.0);
        assert_eq!(m2.determinant(), 0.0);
        assert!(!m2.is_invertible(None));

        let m3 = moebius!(f64, 1.0, 4.0, 0.0, 1.0);
        assert_eq!(m3.determinant(), 1.0);
        assert!(m3.is_invertible(None));

        let m4 = moebius!(i8, 1, 4, 0, -1);
        assert_eq!(m4.determinant(), -1);
        assert!(m4.is_invertible(None));
    }

    #[test]
    fn test_squared_sum() {
        assert_eq!(moebius!(i32, 0, 0, 0, 0).squared_sum(), 0);

        assert_eq!(MoebiusTransformation::<i32>::zero().squared_sum(), 0);

        assert_eq!(
            MoebiusTransformation::<f64>::new(1.0, 0.0, 1.0, 0.0).squared_sum(),
            2.0
        );

        let m1 = moebius!(f32, 1.0, 2.2, 3.0, 4.0);
        assert_eq!((m1 - m1).squared_sum(), 0.0);
    }

    #[test]
    fn test_algebraic_structure() {
        let m = moebius!(f32, 1.0, 2.2, 3.0, 4.0);
        let zero = MoebiusTransformation::<f32>::zero();
        let one = MoebiusTransformation::<f32>::one();
        assert!((m - m).eq(&zero));

        assert!(m + zero == m);
        assert!(zero + m == m);
        assert!(m * zero == zero);
        assert!(zero * m == zero);
        assert!(m * one == m);
        assert!(one * m == m);

        assert!(
            MoebiusTransformation::<i8>::new(1, 1, 0, 0)
                + MoebiusTransformation::<i8>::new(0, 0, 0, 1)
                - MoebiusTransformation::<i8>::new(0, 1, 0, 0)
                == MoebiusTransformation::<i8>::one()
        );

        let m1 = moebius!(f64, -1.0, 2.0, -3.0, 4.0);
        let m2 = moebius!(f64, -5.0, 7.0, 1.0, 5.0);

        assert!(m1 + m2 == MoebiusTransformation::<f64>::new(-6.0, 9.0, -2.0, 9.0));
        assert!(m1 - m2 == MoebiusTransformation::<f64>::new(4.0, -5.0, -4.0, -1.0));
        assert!(m1 * m2 == MoebiusTransformation::<f64>::new(7.0, 3.0, 19.0, -1.0));
    }

    #[test]
    fn test_scalability() {
        assert!(
            MoebiusTransformation::<i8>::new(1, 3, 0, 2) * 3
                == MoebiusTransformation::<i8>::new(3, 9, 0, 6)
        );

        assert!(
            MoebiusTransformation::<f32>::new(2.2, 3.4, -1.2, 0.5) * 2.0
                == MoebiusTransformation::<f32>::new(4.4, 6.8, -2.4, 1.0)
        );
    }

    #[test]
    fn test_preserve_determinant() {
        let m1 = MoebiusTransformation::<i8>::new(1, 4, 0, -1);
        assert_eq!(m1.determinant(), -1);
        assert_eq!((m1 * m1).determinant(), 1);

        let m2 = MoebiusTransformation::<i8>::new(1, 4, -2, -1);
        assert_eq!(m2.determinant(), 7);
        assert_eq!((m2 * m2).determinant(), 49);
        assert_eq!((m1 * m2).determinant(), m1.determinant() * m2.determinant());
        assert_eq!((m1 * m2).determinant(), -7);
    }

    #[test]
    fn test_inverse() {
        let m = MoebiusTransformation::<i8>::new(1, 1, 0, 0);
        assert_eq!(m.determinant(), 0);
        assert!(m.inverse(None).is_none());

        let m = MoebiusTransformation::<f32>::new(1.1, 2.0, -0.5, 3.0);
        assert_eq!(m.determinant(), 4.3);
        assert!(m.inverse(None).is_some());
        assert_eq!(
            m.inverse(None).unwrap().determinant(),
            1.0 / m.determinant()
        );

        assert!(m * m.inverse(None).unwrap() == MoebiusTransformation::<f32>::one());
        assert!(m.inverse(None).unwrap() * m == MoebiusTransformation::<f32>::one());

        // numerical check
        let m = MoebiusTransformation::<f32>::new(
            1.100000000002,
            2.000000000007,
            0.000000000005,
            3.000000000001,
        );
        assert!(m.inverse(None).is_some());
        assert_eq!(
            m.inverse(None).unwrap().determinant(),
            1.0 / m.determinant()
        );

        let numerical_one = m * m.inverse(None).unwrap();
        assert!(numerical_one.is_one(Some(1e-7)));
        assert!(!numerical_one.is_one(Some(1e-8)));
    }
}
