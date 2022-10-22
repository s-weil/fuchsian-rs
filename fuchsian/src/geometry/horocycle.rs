use super::{basics::Mid, boundary::BoundaryPoint};
use crate::{
    algebraic_extensions::{AddIdentity, MulIdentity, Numeric, NumericAddIdentity},
    group_action::{Action, SpecialLinear},
    moebius::MoebiusTransformation,
    NUMERIC_THRESHOLD,
};
use num_complex::Complex;
use std::ops::{Div, Mul};

pub struct HoroCycle<T> {
    /// The `touchpoint` at infinity.
    pub boundary: BoundaryPoint<T>,
    pub height: T,
    pub scale: T,
}

impl<T> HoroCycle<T> {
    pub fn new(height: T) -> Self
    where
        T: MulIdentity,
    {
        Self {
            boundary: BoundaryPoint::Infinity,
            height,
            scale: MulIdentity::one(),
        }
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
            scale: self.scale.clone(),
        }
    }
}

/// Implement Action for Moebius transformations on the boundary.
impl<T> Action<HoroCycle<T>> for MoebiusTransformation<T>
where
    MoebiusTransformation<T>: SpecialLinear<T> + Action<BoundaryPoint<T>>,
    T: Numeric + Copy + NumericAddIdentity + Div<Output = T>,
{
    fn map(&self, x: &HoroCycle<T>) -> HoroCycle<T> {
        match x.boundary {
            BoundaryPoint::Infinity => {
                let scale = if self.c.is_zero(Some(NUMERIC_THRESHOLD)) {
                    x.scale
                } else {
                    x.scale / (self.c * self.c) // / 2.0
                };

                HoroCycle {
                    boundary: self.map(&x.boundary),
                    scale,
                    height: x.height,
                }
            }
            BoundaryPoint::Regular(t) => {
                let denom = self.c * t + self.d;

                if denom.is_zero(Some(NUMERIC_THRESHOLD)) {
                    HoroCycle {
                        boundary: self.map(&x.boundary),
                        scale: x.scale / x.scale, // TODO: require MulIdentity
                        height: x.height,
                    }
                } else {
                    let scale = x.scale / (denom * denom); // TODO: is this correct?
                    HoroCycle {
                        boundary: self.map(&x.boundary),
                        scale,
                        height: x.height,
                    }
                }
            }
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

impl<T> From<HoroCycle<T>> for ParametricHoroCycle<T>
where
    // T: Clone,
    T: Clone + Mul<Output = T> + Div<Output = T> + Mid + AddIdentity,
{
    fn from(hc: HoroCycle<T>) -> Self {
        match hc.boundary {
            BoundaryPoint::Infinity => {
                ParametricHoroCycle::Line(hc.scale.clone() * hc.height.clone())
            }
            BoundaryPoint::Regular(b) => {
                let diameter = hc.scale.clone() / hc.height.clone();
                let zero: T = AddIdentity::zero();
                let radius = zero.mid(&diameter);
                let center = Complex::<T>::new(b, radius.clone());
                ParametricHoroCycle::Circle(EuclieanCircle { center, radius })
            }
        }
    }
}

/* TODO: generalize below


/// The [`HoroCycle`](https://en.wikipedia.org/wiki/Horocycle) in the hyperbolic (Poincare) upper half plane (within C)
/// is either a `Euclidean` circle tangent to the boundary line, i.e. based at a boundary point (the touchpoint), or,
/// the boundary of a half-plane parallel to the real line.
/// `SpecialLinear` preserves horocycles (maps horocycles to horocycles).
/// Given the [`Busemann function`](https://en.wikipedia.org/wiki/Busemann_function) of a boundary point `$\xi$`,
/// the level sets of this function correspond to horocycles based at `$\xi$` and exhaust the hyperbolic space.
/// Conversely, each `horocycle` is the level set of a Busemann function.
///
/// <b>Disclaimer</b>
/// For simplicity, we will use a `height function` in the following which is the Busemann function based at `$\infty$`
/// such that the Horocycle `Im(z) == 1` is the levelset of `0`.
///
pub struct HoroCycle<T> {
    /// The `touchpoint` at infinity.
    pub boundary: BoundaryPoint<T>,
    /// The Busemann distance of the horocycle () SADASDASDASDASD
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
    T: Copy + NumericAddIdentity + Mul<Output = T> + Div<Output = T>,
{
    fn map(&self, x: &HoroCycle<T>) -> HoroCycle<T> {
        let height = if self.c.is_zero(Some(NUMERIC_THRESHOLD)) {
            x.height
        } else {
            // (1.0 - self.c * self.c) / self.c
            // TODO: this is wrong!
            x.height / (self.c * self.c)
        };
        HoroCycle {
            boundary: self.map(&x.boundary),
            height,
        }
    }
}

/ TODO: clarify and correct the relation between height and radius
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

*/
