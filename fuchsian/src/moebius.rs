use std::ops::{Add, Div, Mul, Sub};

/// Mimic some features of a field, without actual structure
pub trait Numeric:
    Sized
    + std::marker::Copy
    + Add<Output = Self>
    + Div<Output = Self>
    + Mul<Output = Self>
    + Sub<Output = Self>
{
    fn is_zero(&self, threshold: Option<Self>) -> bool;
}

#[macro_export]
macro_rules! impl_numeric {
    ($impl_type:ty) => {
        impl Numeric for $impl_type {
            fn is_zero(&self, threshold: Option<Self>) -> bool {
                match threshold {
                    Some(tol) => self.abs() <= tol,
                    None => self.abs() == 0.0,
                }
            }
        }
    };
}

impl_numeric! { f32 }
impl_numeric! { f64 }

/// https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
/// Corresponds to the matrix-vector multiplication of the matrix
/// [[a, b]; [c, d]] * [x, y]
/// for [x, y] in the (real) Euclidean vector space.
pub struct MoebiusTransformation<T>
where
    T: Numeric,
{
    a: T,
    b: T,
    c: T,
    d: T,
}

impl<T> MoebiusTransformation<T>
where
    T: Numeric,
{
    pub fn new(a: T, b: T, c: T, d: T) -> Self {
        Self { a, b, c, d }
    }

    pub fn determinant(&self) -> T {
        self.a * self.d - self.b * self.c
    }

    pub fn is_invertible(&self, threshold: Option<T>) -> bool {
        !self.determinant().is_zero(threshold)
    }
}

#[cfg(test)]
mod tests {
    use super::MoebiusTransformation;

    #[test]
    fn test_is_invertible() {
        let m1 = MoebiusTransformation::<f32>::new(1.0, 2.0, 3.0, 4.0);
        assert_eq!(m1.determinant(), -2.0);
        assert!(m1.is_invertible(None));

        let m2 = MoebiusTransformation::<f64>::new(1.0, 2.0, 2.0, 4.0);
        assert_eq!(m2.determinant(), 0.0);
        assert!(!m2.is_invertible(None));
    }
}
