use super::boundary::BoundaryPoint;
use crate::{
    group_action::{Action, SpecialLinear},
    moebius::MoebiusTransformation,
};

/// The [HoroCycle](https://en.wikipedia.org/wiki/Horocycle) in the hyperbolic (Poincare) upper half plane (within C)
/// is a `Euclidean` circle based at a boundary point (the touchpoint) and determined through its height given
/// in terms of the [Busemann function](https://en.wikipedia.org/wiki/Busemann_function).
pub struct HoroCycle<T> {
    /// The `touchpoint` at infinity.
    pub boundary: BoundaryPoint<T>,
    /// The Busemann distance of the actual horocycle
    pub height: T,
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
