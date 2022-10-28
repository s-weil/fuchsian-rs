use super::{
    basics::{Drawable2d, EuclideanCircle, Mid},
    boundary::BoundaryPoint,
};
use crate::{
    algebraic_extensions::{AddIdentity, MulIdentity, Numeric, NumericAddIdentity},
    group_action::{Action, SpecialLinear},
    moebius::MoebiusTransformation,
};
use num_complex::Complex;
use std::ops::Div;

/// A [`HoroCycle`](https://en.wikipedia.org/wiki/Horocycle) in the hyperbolic space
/// is in general defined as the level set of a [`Busemann function`](https://en.wikipedia.org/wiki/Busemann_function) of a boundary point `$\xi$`,
/// For simplicity we avoid this definition and use of the geometric outcomes in the (Poincare) upper half plane (within C), namely:
/// - $H_t$, a line parallel to the real axis at Euclidean height $t$ corresponding to the level set $Im(z) = t$, or,
/// - an Eucldiean circle tangent to the real axis

pub enum GeometricHorocCycle<T> {
    /// A half-circle with center on the real axis within C
    TangencyCircle(TangencyCircle<T>),
    /// $H_t$ A half-line (level set Im(z) = t) perpendicular to the real axis within C
    Line(T),
}

impl<T> GeometricHorocCycle<T> {
    pub fn new(t: T) -> Self {
        Self::Line(t)
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
}

impl<T> Default for GeometricHorocCycle<T>
where
    T: MulIdentity,
{
    fn default() -> Self {
        GeometricHorocCycle::Line(MulIdentity::one())
    }
}

pub struct TangencyCircle<T> {
    /// The `touchpoint` at infinity.
    boundary: T,
    diameter: T,
}

impl<T> PartialEq for TangencyCircle<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.boundary == other.boundary && self.diameter == other.diameter
    }
}
impl<T> Eq for TangencyCircle<T> where TangencyCircle<T>: PartialEq {}
impl<T> Clone for TangencyCircle<T>
where
    T: Clone,
{
    fn clone(&self) -> TangencyCircle<T> {
        TangencyCircle {
            boundary: self.boundary.clone(),
            diameter: self.diameter.clone(),
        }
    }
}

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

impl<T> PartialEq for GeometricHorocCycle<T>
where
    TangencyCircle<T>: PartialEq,
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (GeometricHorocCycle::Line(t), GeometricHorocCycle::Line(s)) => s == t,
            (GeometricHorocCycle::TangencyCircle(t), GeometricHorocCycle::TangencyCircle(s)) => {
                s == t
            }
            _ => false,
        }
    }
}
impl<T> Eq for GeometricHorocCycle<T> where GeometricHorocCycle<T>: PartialEq {}
impl<T> Clone for GeometricHorocCycle<T>
where
    T: Clone,
    TangencyCircle<T>: Clone,
{
    fn clone(&self) -> GeometricHorocCycle<T> {
        match self {
            GeometricHorocCycle::Line(t) => GeometricHorocCycle::Line(t.clone()),
            GeometricHorocCycle::TangencyCircle(t) => {
                GeometricHorocCycle::TangencyCircle(t.clone())
            }
        }
    }
}

/// Implement Action for Moebius transformations on the boundary.
impl<T> Action<GeometricHorocCycle<T>> for MoebiusTransformation<T>
where
    MoebiusTransformation<T>: SpecialLinear<T> + Action<BoundaryPoint<T>>,
    T: Numeric + Copy + NumericAddIdentity + Div<Output = T> + MulIdentity,
{
    fn map(&self, x: &GeometricHorocCycle<T>) -> GeometricHorocCycle<T> {
        let mapped_boundary = self.map(&x.boundary_point());
        match x {
            GeometricHorocCycle::Line(t) => {
                // See Dal'bo Lemma 3.19, page 30
                match mapped_boundary {
                    BoundaryPoint::Infinity => GeometricHorocCycle::Line(*t), // case c=0, must be a horocyclic isometry already
                    BoundaryPoint::Regular(boundary) => {
                        let one: T = MulIdentity::one();
                        let diameter = one / (*t * self.c * self.c);
                        let circle = TangencyCircle { boundary, diameter };
                        GeometricHorocCycle::TangencyCircle(circle)
                    }
                }
            }
            GeometricHorocCycle::TangencyCircle(circle) => {
                match mapped_boundary {
                    BoundaryPoint::Infinity => {
                        let h = circle.diameter * self.c * self.c;
                        let one: T = MulIdentity::one();
                        GeometricHorocCycle::Line(one / h)
                    }
                    BoundaryPoint::Regular(boundary) => {
                        // case denom != 0
                        let denom = self.c * circle.boundary + self.d;
                        let diameter = circle.diameter / (denom * denom);
                        let circle = TangencyCircle { boundary, diameter };
                        GeometricHorocCycle::TangencyCircle(circle)
                    }
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
