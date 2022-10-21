use num_complex::Complex;

use super::boundary::BoundaryPoint;
use crate::{
    group_action::{Action, SpecialLinear},
    moebius::MoebiusTransformation,
};

/// The [HoroCycle](https://en.wikipedia.org/wiki/Horocycle) in the hyperbolic (Poincare) upper half plane (within C)
/// is either a `Euclidean` circle based at a boundary point (the touchpoint) and determined through its height given
/// in terms of the [Busemann function](https://en.wikipedia.org/wiki/Busemann_function), or,
/// the boundary of a half-plane parallel to the real line.
/// SpecialLinear preserves horocycles (maps horocycles to horocycles).
pub struct HoroCycle<T> {
    /// The `touchpoint` at infinity.
    pub boundary: BoundaryPoint<T>,
    /// The Busemann distance of the actual horocycle
    pub height: T,
}

impl<T> HoroCycle<T> {
    pub fn new(boundary: BoundaryPoint<T>, height: T) -> Self {
        Self { boundary, height }
    }
}

impl<T> PartialEq for HoroCycle<T>
where
    BoundaryPoint<T>: PartialEq,
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.boundary == other.boundary && self.height == other.height
    }
}
impl<T> Eq for HoroCycle<T>
where
    BoundaryPoint<T>: PartialEq,
    T: PartialEq,
{
}
impl<T> Copy for HoroCycle<T> where T: Copy {}
impl<T> Clone for HoroCycle<T>
where
    T: Clone,
{
    fn clone(&self) -> HoroCycle<T> {
        HoroCycle {
            boundary: self.boundary.clone(),
            height: self.height.clone(),
        }
    }
}

/// Implement Action for Moebius transformations on the boundary.
impl<T> Action<HoroCycle<T>> for MoebiusTransformation<T>
where
    MoebiusTransformation<T>: SpecialLinear<T> + Action<BoundaryPoint<T>>,
    T: Copy,
{
    fn map(&self, x: &HoroCycle<T>) -> HoroCycle<T> {
        HoroCycle {
            boundary: self.map(&x.boundary),
            height: x.height,
        }
    }
}

pub enum ParametricHoroCycle<T> {
    /// A half-circle with center on the real axis within C
    Circle(EuclieanCircle<T>),
    /// A half-line perpendicular to the real axis within C
    Line(T), // height / distance to the real line
}

/// The parametrization of a Euclidean circle in C.
pub struct EuclieanCircle<T> {
    pub center: Complex<T>,
    pub radius: T,
}

impl<T> PartialEq for EuclieanCircle<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.center == other.center && self.radius == other.radius
    }
}
impl<T> Eq for EuclieanCircle<T> where T: PartialEq {}

// TODO: clarify and correct the relation between height and radius
impl<T> From<HoroCycle<T>> for ParametricHoroCycle<T>
where
    T: Clone,
{
    fn from(hc: HoroCycle<T>) -> Self {
        match hc.boundary {
            BoundaryPoint::Infinity => ParametricHoroCycle::Line(hc.height),
            BoundaryPoint::Regular(b) => {
                let radius = hc.height.clone(); // sqrt / 2?
                let center = Complex::<T>::new(b, radius.clone());
                ParametricHoroCycle::Circle(EuclieanCircle { center, radius })
            }
        }
    }
}
