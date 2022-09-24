use crate::algebraic_extensions::{AddIdentity, MulIdentity, Numeric, NumericAddIdentity};
use std::{
    fmt,
    ops::{Add, Mul, Neg, Sub},
};

/// https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
/// Corresponds to the matrix-vector multiplication of the matrix
/// [[a, b]; [c, d]] * [x, y]
/// for [x, y] in the (real) Euclidean vector space.
pub struct MoebiusTransformation<T> {
    a: T,
    b: T,
    c: T,
    d: T,
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
        self.a * self.d - self.b * self.c
    }

    pub fn is_invertible(&self, threshold: Option<f64>) -> bool
    where
        T: Numeric + NumericAddIdentity + std::marker::Copy,
    {
        !self.determinant().is_zero(threshold)
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

// ########################
// Establish algebraic structure for Moebius-transformations
// ########################

impl<T> AddIdentity for MoebiusTransformation<T>
where
    T: AddIdentity,
{
    fn zero() -> Self {
        MoebiusTransformation::new(T::zero(), T::zero(), T::zero(), T::zero())
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

// ########################
// Establish vector space structure for Moebius-transformations
// ########################

// addition see above

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

#[cfg(test)]
mod tests {
    use super::MoebiusTransformation;
    use crate::algebraic_extensions::{AddIdentity, MulIdentity};

    #[test]
    fn test_is_invertible() {
        let m1 = MoebiusTransformation::<f32>::new(1.0, 2.0, 3.0, 4.0);
        assert_eq!(m1.determinant(), -2.0);
        assert!(m1.is_invertible(None));

        let m2 = MoebiusTransformation::<f64>::new(1.0, 2.0, 2.0, 4.0);
        assert_eq!(m2.determinant(), 0.0);
        assert!(!m2.is_invertible(None));

        let m3 = MoebiusTransformation::<f64>::new(1.0, 4.0, 0.0, 1.0);
        assert_eq!(m3.determinant(), 1.0);
        assert!(m3.is_invertible(None));

        let m4 = MoebiusTransformation::<i8>::new(1, 4, 0, -1);
        assert_eq!(m4.determinant(), -1);
        assert!(m4.is_invertible(None));
    }

    #[test]
    fn test_squared_sum() {
        assert_eq!(
            MoebiusTransformation::<i32>::new(0, 0, 0, 0).squared_sum(),
            0
        );

        assert_eq!(MoebiusTransformation::<i32>::zero().squared_sum(), 0);

        assert_eq!(
            MoebiusTransformation::<f64>::new(1.0, 0.0, 1.0, 0.0).squared_sum(),
            2.0
        );

        let m1 = MoebiusTransformation::<f32>::new(1.0, 2.2, 3.0, 4.0);
        assert_eq!((m1 - m1).squared_sum(), 0.0);
    }

    #[test]
    fn test_group_structure() {
        let m = MoebiusTransformation::<f32>::new(1.0, 2.2, 3.0, 4.0);
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

        let m1 = MoebiusTransformation::<f64>::new(-1.0, 2.0, -3.0, 4.0);
        let m2 = MoebiusTransformation::<f64>::new(-5.0, 7.0, 1.0, 5.0);

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
}
