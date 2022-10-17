use super::{
    basics::{Distance, Mid},
    boundary::BoundaryPoint,
};
use crate::{
    group_action::{Action, SpecialLinear},
    moebius::MoebiusTransformation,
    set_extensions::SetRestriction,
};


/// An oriented geodesic (from `start` to `end`) in the hyperbolic upper half plane (within C) is uniquely determined by
/// its <i>two distinct</i> endpoints on the boundary. Geometrically it will be either an Euclidean (half) line perpendicular to the real axis,
/// or, an Euclidean arc (half circle) with center on the real axis.
/// SpecialLinear preserves geodescis (maps geodesics to geodesics).
pub struct GeodesicBoundary<T> {
    start: BoundaryPoint<T>,
    end: BoundaryPoint<T>,
}

impl<T> SetRestriction for GeodesicBoundary<T>
where
    T: PartialEq,
{
    fn condition(&self) -> bool {
        self.start != self.end
    }
}

impl<T> GeodesicBoundary<T>
where
    T: PartialEq,
{
    pub fn end_points(&self) -> (BoundaryPoint<T>, BoundaryPoint<T>)
    where
        T: Clone,
    {
        (self.start.clone(), self.end.clone())
    }
}

/// Implement `Action` for Moebius transformations on Geodesics.
impl<T> Action<GeodesicBoundary<T>> for MoebiusTransformation<T>
where
    MoebiusTransformation<T>: SpecialLinear<T> + Action<BoundaryPoint<T>>,
{
    fn map(&self, x: &GeodesicBoundary<T>) -> GeodesicBoundary<T> {
        GeodesicBoundary {
            start: self.map(&x.start),
            end: self.map(&x.end),
        }
    }
}

/// Geodesics of the upper half plane (within C).
/// SpecialLinear preserves the boundary (maps geodesics to geodesics).
pub enum GeodesicLine<T> {
    /// A half-circle with center on the real axis within C
    Arc(Arc<T>),
    /// A half-line perpendicular to the real axis within C
    Line(T), // touchpoint
}

/// The parametrization of a geodesic line in case of an arc.
pub struct Arc<T> {
    center: T,
    radius: T,
}

impl<T> Arc<T> {
    pub fn from_antipodal(p1: T, p2: T) -> Self
    where
        T: Distance<T> + Mid,
    {
        // let distance = p1.dist(&p2);
        // let radius = T::zero().mid(distance);
        let center = p1.mid(&p2);
        let radius = p1.dist(&center);
        Self { center, radius }
    }
}

// /// NOTE: for T in { i8, i32, i64 } etc, there is in general NO unique geodesic parametrization
impl<T> From<GeodesicBoundary<T>> for GeodesicLine<T>
where
    T: Distance<T> + Mid,
{
    fn from(gb: GeodesicBoundary<T>) -> Self {
        // if &gb.start == &gb.end {
        //     panic!("not possible");
        // }
        match (gb.start, gb.end) {
            (BoundaryPoint::Infinity, BoundaryPoint::Regular(t))
            | (BoundaryPoint::Regular(t), BoundaryPoint::Infinity) => GeodesicLine::Line(t),
            (BoundaryPoint::Regular(s), BoundaryPoint::Regular(e)) => {
                let arc = Arc::from_antipodal(s, e);
                GeodesicLine::Arc(arc)
            }
            _ => panic!("not possible"),
        }
    }
}
