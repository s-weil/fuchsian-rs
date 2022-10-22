use super::{
    basics::{Distance, Mid},
    boundary::BoundaryPoint,
};
use crate::{
    group_action::{Action, SpecialLinear},
    moebius::MoebiusTransformation,
    set_extensions::SetRestriction,
};

/// An `oriented geodesic` (from `start` to `end`) in the hyperbolic upper half plane (within C) is uniquely determined by
/// its <i>two distinct</i> endpoints on the boundary. Geometrically it will be either an Euclidean (half) line perpendicular to the real axis,
/// or, an Euclidean arc (half circle) with center on the real axis.
/// `SpecialLinear` preserves geodescis (maps geodesics to geodesics).
pub struct GeodesicBoundary<T> {
    start: BoundaryPoint<T>,
    end: BoundaryPoint<T>,
}

impl<T> GeodesicBoundary<T> {
    pub fn new(start: BoundaryPoint<T>, end: BoundaryPoint<T>) -> Self {
        Self { start, end }
    }
}

impl<T> SetRestriction for GeodesicBoundary<T>
where
    T: PartialEq,
{
    fn condition(&self) -> bool {
        self.start != self.end
    }
}

impl<T> PartialEq for GeodesicBoundary<T>
where
    BoundaryPoint<T>: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.start == other.end && self.end == other.end
    }
}
impl<T> Eq for GeodesicBoundary<T> where BoundaryPoint<T>: PartialEq {}
impl<T> Copy for GeodesicBoundary<T> where T: Copy {}
impl<T> Clone for GeodesicBoundary<T>
where
    T: Clone,
{
    fn clone(&self) -> GeodesicBoundary<T> {
        GeodesicBoundary {
            start: self.start.clone(),
            end: self.end.clone(),
        }
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

impl<T> PartialEq for GeodesicLine<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (GeodesicLine::Line(t1), GeodesicLine::Line(t2)) => t1 == t2,
            (GeodesicLine::Arc(a1), GeodesicLine::Arc(a2)) => a1 == a2,
            _ => false,
        }
    }
}

/// The parametrization of a geodesic line in case of an arc.
pub struct Arc<T> {
    pub center: T,
    pub radius: T,
}

impl<T> PartialEq for Arc<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.center == other.center && self.radius == other.radius
    }
}
impl<T> Eq for Arc<T> where T: PartialEq {}

impl<T> Arc<T> {
    pub fn from_antipodal(p1: T, p2: T) -> Self
    where
        T: Distance<T> + Mid,
    {
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

#[cfg(test)]
mod tests {
    use super::{GeodesicBoundary, GeodesicLine};
    use crate::{
        fuchsian_group::{FuchsianGroup, SpecialLinearMoebiusTransformation},
        geometry::{boundary::BoundaryPoint, geodesics::Arc},
        group_action::{Action, Orbit},
        moebius,
        moebius::MoebiusTransformation,
        set_extensions::SetRestriction,
    };

    #[test]
    fn create() {
        let g1 = GeodesicBoundary {
            start: BoundaryPoint::Regular(1.0),
            end: BoundaryPoint::Infinity,
        };
        let l1 = GeodesicLine::Line(1.0_f64);
        assert!(GeodesicLine::from(g1) == l1);

        let g2 = GeodesicBoundary {
            start: BoundaryPoint::Regular(2.0),
            end: BoundaryPoint::Regular(-2.0),
        };
        let l2 = GeodesicLine::Arc(Arc {
            center: 0.0,
            radius: 2.0,
        });
        assert!(GeodesicLine::from(g2) == l2);
    }

    #[test]
    fn hyperbolic_map_geodesic() {
        let sl =
            SpecialLinearMoebiusTransformation::try_from(moebius!(f64, 2.0, 0.0, 0.0, 0.5), None)
                .unwrap();

        let g = GeodesicBoundary {
            start: BoundaryPoint::Regular(1.0),
            end: BoundaryPoint::Regular(-1.0),
        };
        let mg = sl.map(&g);
        let l = GeodesicLine::Arc(Arc {
            center: 0.0,
            radius: 4.0,
        });
        assert!(GeodesicLine::from(mg) == l);

        let g = GeodesicBoundary {
            start: BoundaryPoint::Infinity,
            end: BoundaryPoint::Regular(-1.0),
        };
        let mg = sl.map(&g);
        let l = GeodesicLine::Line(-4.0);
        assert!(GeodesicLine::from(mg) == l);
    }

    #[test]
    fn parabolic_map_geodesic() {
        let sl =
            SpecialLinearMoebiusTransformation::try_from(moebius!(f64, 1.0, 3.0, 0.0, 1.0), None)
                .unwrap();

        let g = GeodesicBoundary {
            start: BoundaryPoint::Regular(1.0),
            end: BoundaryPoint::Regular(-1.0),
        };
        let mg = sl.map(&g);
        let l = GeodesicLine::Arc(Arc {
            center: 3.0,
            radius: 1.0,
        });
        assert!(GeodesicLine::from(mg) == l);

        let g = GeodesicBoundary {
            start: BoundaryPoint::Infinity,
            end: BoundaryPoint::Regular(-1.0),
        };
        let mg = sl.map(&g);
        let l = GeodesicLine::Line(2.0);
        assert!(GeodesicLine::from(mg) == l);
    }

    #[test]
    fn test_geodesic_orbit_modular_group() {
        // see https://en.wikipedia.org/wiki/Modular_group
        // the modular group is generated by the transformations z -> z+1 and z -> -1/z,
        // corresponding to the Moebius transformations g = [ 1, 1; 0, 1 ] and h = [ 0, -1; 1, 0 ]
        let g = MoebiusTransformation::<f32>::new(1.0, 1.0, 0.0, 1.0);
        let h = MoebiusTransformation::<f32>::new(0.0, -1.0, 1.0, 0.0);

        let fuchsian_group = FuchsianGroup::create_from_valid(vec![g, h]);

        let geodesic_arc = GeodesicBoundary::try_new(GeodesicBoundary {
            start: BoundaryPoint::Regular(-1.0),
            end: BoundaryPoint::Regular(1.0),
        })
        .unwrap();

        let orbit = Orbit::sample(&fuchsian_group, &geodesic_arc, 100, None);

        assert_eq!(orbit.points.len(), 100);
    }
}
