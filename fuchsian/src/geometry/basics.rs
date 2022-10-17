use crate::set_extensions::SetRestriction;

pub trait UpperHalfSpace: SetRestriction {
    fn condition(&self) -> bool;
}

pub trait Distance<T> {
    fn dist(&self, other: &Self) -> T;
}

impl Distance<f32> for f32 {
    fn dist(&self, other: &Self) -> f32 {
        (self - other).abs()
    }
}

impl Distance<f64> for f64 {
    fn dist(&self, other: &Self) -> f64 {
        (self - other).abs()
    }
}

// TODO: impl dist for complex, hyperbolic, boundary points

/// NOTE: for T in { i8, i32, i64 } etc, there is in general NO unique mid point
pub trait Mid {
    fn mid(&self, other: &Self) -> Self;
}

impl Mid for f32 {
    fn mid(&self, other: &Self) -> Self {
        (self + other) / 2.0
    }
}

impl Mid for f64 {
    fn mid(&self, other: &Self) -> Self {
        (self + other) / 2.0
    }
}

// TODO: impl Mid for complex, hyperbolic points
