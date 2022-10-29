use super::{
    basics::{Drawable2d, EuclideanCircle, Mid},
    boundary::BoundaryPoint,
};
use crate::{
    algebraic_extensions::{AddIdentity, MulIdentity, Numeric, NumericAddIdentity},
    group_action::{Action, SpecialLinear},
    moebius::MoebiusTransformation,
    NUMERIC_THRESHOLD,
};
use num_complex::Complex;
use std::ops::Div;

/// A [`HoroCycle`](https://en.wikipedia.org/wiki/Horocycle) in the hyperbolic space
/// is in general defined as the level set of a [`Busemann function`](https://en.wikipedia.org/wiki/Busemann_function) of a boundary point `$\xi$`,
/// For simplicity we avoid this definition and use of the geometric outcomes in the (Poincare) upper half plane (within C), namely:
/// - $H_t$, a line parallel to the real axis at Euclidean height $t$ corresponding to the level set $Im(z) = t$, or,
/// - an Eucldiean circle tangent to the real axis
#[derive(Debug, PartialEq, Eq, Clone)]
pub enum GeometricHorocCycle<T> {
    /// A half-circle with center on the real axis within C
    TangencyCircle(TangencyCircle<T>),
    /// $H_t$ A half-line (level set Im(z) = t) perpendicular to the real axis within C
    Line(T),
}

impl<T> GeometricHorocCycle<T> {
    pub fn new(boundary: BoundaryPoint<T>, height_or_diam: T) -> Self {
        // if height_or_diam <= 0.0 {
        //     panic!("Provided non-positive height or diameter");
        // }
        match boundary {
            BoundaryPoint::Infinity => GeometricHorocCycle::Line(height_or_diam),
            BoundaryPoint::Regular(t) => GeometricHorocCycle::TangencyCircle(TangencyCircle {
                boundary: t,
                diameter: height_or_diam,
            }),
        }
    }

    pub fn new_line(t: T) -> Self {
        Self::new(BoundaryPoint::Infinity, t)
    }

    pub fn boundary_point(&self) -> BoundaryPoint<T>
    where
        T: Clone,
    {
        match self {
            Self::Line(_) => BoundaryPoint::Infinity,
            Self::TangencyCircle(TangencyCircle { boundary, .. }) => {
                BoundaryPoint::Regular(boundary.clone())
            }
        }
    }

    pub fn height_or_diameter(&self) -> T
    where
        T: Copy,
    {
        match self {
            GeometricHorocCycle::Line(t) => *t,
            GeometricHorocCycle::TangencyCircle(TangencyCircle { diameter, .. }) => *diameter,
        }
    }
}

impl<T> Default for GeometricHorocCycle<T>
where
    T: MulIdentity,
{
    fn default() -> Self {
        GeometricHorocCycle::Line(MulIdentity::one())
    }
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct TangencyCircle<T> {
    /// The `touchpoint` at infinity.
    boundary: T,
    diameter: T,
}

// impl<T> PartialEq for TangencyCircle<T>
// where
//     T: PartialEq,
// {
//     fn eq(&self, other: &Self) -> bool {
//         self.boundary == other.boundary && self.diameter == other.diameter
//     }
// }
// impl<T> Eq for TangencyCircle<T> where TangencyCircle<T>: PartialEq {}
// impl<T> Clone for TangencyCircle<T>
// where
//     T: Clone,
// {
//     fn clone(&self) -> TangencyCircle<T> {
//         TangencyCircle {
//             boundary: self.boundary.clone(),
//             diameter: self.diameter.clone(),
//         }
//     }
// }

impl<T> From<&TangencyCircle<T>> for EuclideanCircle<T>
where
    T: Clone + Mid + AddIdentity,
{
    fn from(circle: &TangencyCircle<T>) -> Self {
        let zero: T = AddIdentity::zero();
        let radius = zero.mid(&circle.diameter);
        let center = Complex::<T>::new(circle.boundary.clone(), radius.clone());
        EuclideanCircle { center, radius }
    }
}

// impl<T> PartialEq for GeometricHorocCycle<T>
// where
//     TangencyCircle<T>: PartialEq,
//     T: PartialEq,
// {
//     fn eq(&self, other: &Self) -> bool {
//         match (self, other) {
//             (GeometricHorocCycle::Line(t), GeometricHorocCycle::Line(s)) => s == t,
//             (GeometricHorocCycle::TangencyCircle(t), GeometricHorocCycle::TangencyCircle(s)) => {
//                 s == t
//             }
//             _ => false,
//         }
//     }
// }
// impl<T> Eq for GeometricHorocCycle<T> where GeometricHorocCycle<T>: PartialEq {}
// impl<T> Clone for GeometricHorocCycle<T>
// where
//     T: Clone,
//     TangencyCircle<T>: Clone,
// {
//     fn clone(&self) -> GeometricHorocCycle<T> {
//         match self {
//             GeometricHorocCycle::Line(t) => GeometricHorocCycle::Line(t.clone()),
//             GeometricHorocCycle::TangencyCircle(t) => {
//                 GeometricHorocCycle::TangencyCircle(t.clone())
//             }
//         }
//     }
// }

/// Implement Action for Moebius transformations on the boundary.
impl<T> Action<GeometricHorocCycle<T>> for MoebiusTransformation<T>
where
    MoebiusTransformation<T>: SpecialLinear<T> + Action<BoundaryPoint<T>>,
    T: Numeric + Copy + NumericAddIdentity + Div<Output = T> + MulIdentity,
{
    fn map(&self, x: &GeometricHorocCycle<T>) -> GeometricHorocCycle<T> {
        match x {
            GeometricHorocCycle::Line(t) => {
                // See Dal'bo Lemma 3.19, page 30
                if self.c.is_zero(Some(NUMERIC_THRESHOLD)) {
                    // BoundaryPoint::infty is one fixed point
                    if (self.a + (-self.d)).is_zero(Some(NUMERIC_THRESHOLD)) {
                        GeometricHorocCycle::Line(*t)
                    } else {
                        // let other_fixed_point = self.a * self.b  / (1 - self.a * self.a);
                        let height = self.a * self.a * *t;
                        GeometricHorocCycle::Line(height)
                    }
                } else {
                    let one: T = MulIdentity::one();
                    let diameter = one / (*t * self.c * self.c);
                    let circle = TangencyCircle {
                        boundary: self.a / self.c,
                        diameter,
                    };
                    GeometricHorocCycle::TangencyCircle(circle)
                }
            }
            GeometricHorocCycle::TangencyCircle(circle) => {
                let denom = self.c * circle.boundary + self.d;
                if denom.is_zero(Some(NUMERIC_THRESHOLD)) {
                    let h = circle.diameter * self.c * self.c;
                    let one: T = MulIdentity::one();
                    GeometricHorocCycle::Line(one / h)
                } else {
                    let diameter = circle.diameter / (denom * denom); // TODO: check if correct
                    let boundary = (self.a * circle.boundary + self.b) / denom;
                    let circle = TangencyCircle { boundary, diameter };
                    GeometricHorocCycle::TangencyCircle(circle)
                }
            }
        }
    }
}

impl Drawable2d<f64> for GeometricHorocCycle<f64> {
    fn draw(&self, n_curve_points: usize) -> Vec<(f64, f64)> {
        match self {
            GeometricHorocCycle::Line(height) => {
                vec![
                    (-(n_curve_points as f64), *height),
                    (n_curve_points as f64, *height),
                ]
            }
            GeometricHorocCycle::TangencyCircle(tangency_circle) => {
                let eucl_circle = EuclideanCircle::from(tangency_circle);
                eucl_circle.draw(n_curve_points)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        geometry::{boundary::BoundaryPoint, horocycle::GeometricHorocCycle},
        group_action::Action,
        moebius::MoebiusTransformation,
    };
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_action_horocyclic() {
        let h = MoebiusTransformation::<f64>::new(1.0, 10.0, 0.0, 1.0);
        let dflt_horocycle = GeometricHorocCycle::default();
        let mapped_hc = h.map(&dflt_horocycle);
        assert_eq!(mapped_hc, dflt_horocycle);

        let h_inv = h.inverse(None).unwrap();
        let remapped_hc = h_inv.map(&mapped_hc);
        assert_eq!(remapped_hc, dflt_horocycle);

        let tang_hc = GeometricHorocCycle::new(BoundaryPoint::Regular(0.0), 1.0);
        let mapped_hc = h.map(&tang_hc);
        assert_eq!(
            mapped_hc,
            GeometricHorocCycle::new(BoundaryPoint::Regular(10.0), 1.0)
        );
    }

    #[test]
    fn test_action_elliptic() {
        let h = MoebiusTransformation::<f64>::new(0.0, -1.0, 1.0, 0.0);
        let dflt_horocycle = GeometricHorocCycle::default();
        let mapped_hc = h.map(&dflt_horocycle);
        assert_ne!(mapped_hc, dflt_horocycle);
        assert_eq!(mapped_hc.boundary_point(), BoundaryPoint::Regular(0.0));

        let h_inv = h.inverse(None).unwrap();
        let remapped_hc = h_inv.map(&mapped_hc);
        assert_eq!(remapped_hc, dflt_horocycle);

        let tang_hc = GeometricHorocCycle::new(BoundaryPoint::Regular(0.0), 1.0);
        let mapped_hc = h.map(&tang_hc);
        assert_eq!(mapped_hc, GeometricHorocCycle::default());
    }

    #[test]
    fn test_action_hyperbolic() {
        let h = MoebiusTransformation::<f64>::new(5.0, 0.0, 0.0, 0.2);
        let h_inv = h.inverse(None).unwrap();

        let dflt_horocycle = GeometricHorocCycle::default();
        let mapped_hc = h.map(&dflt_horocycle);
        assert_eq!(mapped_hc, GeometricHorocCycle::Line(25.0));

        let remapped_hc = h_inv.map(&mapped_hc);
        assert_eq!(
            remapped_hc.boundary_point(),
            dflt_horocycle.boundary_point()
        );
        assert_abs_diff_eq!(
            remapped_hc.height_or_diameter(),
            dflt_horocycle.height_or_diameter(),
            epsilon = f64::EPSILON
        );

        let tang_hc = GeometricHorocCycle::new(BoundaryPoint::Regular(0.0), 1.0);
        let mapped_hc = h.map(&tang_hc);
        assert_eq!(mapped_hc.boundary_point(), BoundaryPoint::Regular(0.0));
        assert_abs_diff_eq!(mapped_hc.height_or_diameter(), 25.0, epsilon = 1e-12);

        let remapped_hc = h_inv.map(&mapped_hc);
        assert_eq!(remapped_hc.boundary_point(), BoundaryPoint::Regular(0.0));
        assert_abs_diff_eq!(remapped_hc.height_or_diameter(), 1.0, epsilon = 1e-12);

        let tang_hc = GeometricHorocCycle::new(BoundaryPoint::Regular(1.0), 1.0);
        let mapped_hc = h.map(&tang_hc);
        assert_eq!(mapped_hc.boundary_point(), BoundaryPoint::Regular(25.0));
        assert_ne!(mapped_hc.height_or_diameter(), 25.0);
        // TODO: what is the diameter?
        let remapped_hc = h_inv.map(&mapped_hc);
        assert_eq!(remapped_hc.boundary_point(), BoundaryPoint::Regular(1.0));
        assert_abs_diff_eq!(remapped_hc.height_or_diameter(), 1.0, epsilon = 1e-12);
    }
}
