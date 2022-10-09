use super::boundary::BoundaryPoint;
use crate::{
    algebraic_extensions::{Numeric, NumericAddIdentity},
    group_action::{Action, SpecialLinear},
    moebius::MoebiusTransformation,
};
use std::ops::Div;

/// An oriented geodesic (from `start` to `end`) in the hyperbolic upper half plane (within C) is uniquely determined by
/// its <i>two distinct</i> endpoints on the boundary.
/// SpecialLinear preserves geodescis (maps geodesics to geodesics).
pub struct GeodesicBoundary<T> {
    start: BoundaryPoint<T>,
    end: BoundaryPoint<T>,
}

impl<T> GeodesicBoundary<T>
where
    T: PartialEq,
{
    pub fn new(start: BoundaryPoint<T>, end: BoundaryPoint<T>) -> Option<Self> {
        if start != end {
            return Some(Self { start, end });
        }
        None
    }
}

// /// Implement `Action` for Moebius transformations on Geodesics.
// impl<T> Action<GeodesicBoundary<T>> for MoebiusTransformation<T>
// where
//     T: Numeric + Div<Output = T> + NumericAddIdentity + Copy,
//     MoebiusTransformation<T>: SpecialLinear<T>,
// {
//     fn map(&self, x: &GeodesicBoundary<T>) -> GeodesicBoundary<T> {
//         GeodesicBoundary {
//             start: self.map(&x.start),
//             end: self.map(&x.end),
//         }
//     }
// }

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

pub struct Arc<T> {
    center: T,
    radius: T,
}

// impl<T> From<GeodesicBoundary<T>> for GeodesicLine<T>
// where
//     T: Sub<Output = T> + PartialEq,
// {
//     fn from(gb: GeodesicBoundary<T>) -> Self {
//         if &gb.start == &gb.end {
//             panic!("not possible");
//         }
//         match (gb.start, gb.end) {
//             (BoundaryPoint::Infinity, BoundaryPoint::Regular(t))
//             | (BoundaryPoint::Regular(t), BoundaryPoint::Infinity) => GeodesicLine::Line(t),
//             (BoundaryPoint::Regular(s), BoundaryPoint::Regular(e)) => {
//                 let radius = (s - e).abs() / 2.0;
//                 let center = (s.min(e)) + radius;
//                 let arc = Arc { center, radius };
//                 GeodesicLine::Arc(arc)
//             }
//             _ => panic!("not possible"),
//         }
//     }
// }
