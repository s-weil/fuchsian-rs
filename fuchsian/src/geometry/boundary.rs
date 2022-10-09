use std::ops::Div;

use crate::{
    algebraic_extensions::{Numeric, NumericAddIdentity},
    group_action::Action,
    moebius::MoebiusTransformation,
};

pub enum BoundaryPoint<T> {
    Infinity,
    Regular(T),
}

// TODO: action of moebius

const NUMERIC_THRESHOLD: f64 = 1e-16;

/// Implement Action for Moebius transformations on the boundary.
impl<T> Action<BoundaryPoint<T>> for MoebiusTransformation<T>
where
    T: Numeric + Div<Output = T> + NumericAddIdentity + Copy + PartialEq,
{
    fn map(&self, x: &BoundaryPoint<T>) -> BoundaryPoint<T> {
        match x {
            BoundaryPoint::Infinity => {
                if self.c.is_zero(Some(NUMERIC_THRESHOLD)) {
                    BoundaryPoint::Infinity
                } else {
                    BoundaryPoint::Regular(self.a / self.c)
                }
            }
            BoundaryPoint::Regular(t) => {
                let denom = self.c * *t + self.d;
                if denom.is_zero(Some(NUMERIC_THRESHOLD)) {
                    BoundaryPoint::Infinity
                } else {
                    BoundaryPoint::Regular((self.a * *t + self.b) / denom)
                }
            }
        }
    }
}
