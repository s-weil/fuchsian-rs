use crate::{
    algebraic_extensions::{Numeric, NumericAddIdentity},
    group_action::{Action, SpecialLinear},
    moebius::MoebiusTransformation,
    NUMERIC_THRESHOLD,
};
use std::{
    fmt::{self, Debug, Display},
    ops::Div,
};

/// Boundary points of the hyperbolic (Poincare) upper half plane (within C).
/// SpecialLinear preserves the boundary (maps a boundary points to the boundary).
pub enum BoundaryPoint<T> {
    Infinity,
    Regular(T),
}

impl<T> PartialEq for BoundaryPoint<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (BoundaryPoint::Infinity, BoundaryPoint::Infinity) => true,
            (BoundaryPoint::Regular(t1), BoundaryPoint::Regular(t2)) => t1 == t2,
            _ => false,
        }
    }
}
impl<T> Eq for BoundaryPoint<T> where T: PartialEq {}
impl<T> Display for BoundaryPoint<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BoundaryPoint::Infinity => write!(f, "∞"),
            BoundaryPoint::Regular(t) => write!(f, "{}", t),
        }
    }
}
impl<T> Debug for BoundaryPoint<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BoundaryPoint::Infinity => write!(f, "∞"),
            BoundaryPoint::Regular(t) => write!(f, "δ({})", t),
        }
    }
}

/// Implement Action for Moebius transformations on the boundary.
impl<T> Action<BoundaryPoint<T>> for MoebiusTransformation<T>
where
    T: Numeric + Div<Output = T> + NumericAddIdentity + Copy,
    MoebiusTransformation<T>: SpecialLinear<T>,
{
    fn map(&self, x: &BoundaryPoint<T>) -> BoundaryPoint<T> {
        match x {
            BoundaryPoint::Infinity => {
                // Note: a and c cannot be 0 at the same time
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

#[cfg(test)]
mod tests {
    use super::BoundaryPoint;
    use crate::{group_action::Action, moebius::MoebiusTransformation};

    #[test]
    fn test_action_horocyclic() {
        let h = MoebiusTransformation::<f64>::new(1.0, 10.0, 0.0, 1.0);

        let boundary_infty: BoundaryPoint<f64> = BoundaryPoint::Infinity;
        let boundary_regular = BoundaryPoint::Regular(0.0);

        assert_eq!(h.map(&boundary_infty), BoundaryPoint::Infinity);
        assert_eq!(h.map(&boundary_regular), BoundaryPoint::Regular(10.0));
    }

    #[test]
    fn test_action_hyperbolic() {
        let h = MoebiusTransformation::<f64>::new(5.0, 0.0, 0.0, 0.2);

        // fixed point 1
        let boundary_infty: BoundaryPoint<f64> = BoundaryPoint::Infinity;
        assert_eq!(h.map(&boundary_infty), BoundaryPoint::Infinity);

        // fixed point 2
        let boundary_regular = BoundaryPoint::Regular(0.0);
        assert_eq!(h.map(&boundary_regular), BoundaryPoint::Regular(0.0));

        // no fixed point
        let boundary_regular = BoundaryPoint::Regular(1.0);
        assert_eq!(h.map(&boundary_regular), BoundaryPoint::Regular(25.0));
        let boundary_regular = BoundaryPoint::Regular(-1.0);
        assert_eq!(h.map(&boundary_regular), BoundaryPoint::Regular(-25.0));
    }

    // TODO: incorrect
    #[test]
    fn test_action_rotation() {
        // no fixed points
        // rotation by pi/2 = 90 [cos(.), -sin(.); sin(.), cos(.)]
        let h = MoebiusTransformation::<f64>::new(0.0, -1.0, 1.0, 0.0);

        let b = BoundaryPoint::Regular(0.0);
        assert_eq!(h.map(&b), BoundaryPoint::Infinity);

        let b = BoundaryPoint::Infinity;
        assert_eq!(h.map(&b), BoundaryPoint::Regular(0.0));

        let b = BoundaryPoint::Regular(1.0);
        assert_eq!(h.map(&b), BoundaryPoint::Regular(-1.0));

        let b = BoundaryPoint::Regular(-1.0);
        assert_eq!(h.map(&b), BoundaryPoint::Regular(1.0));
    }
}
