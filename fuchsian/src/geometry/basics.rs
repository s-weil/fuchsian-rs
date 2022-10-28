use num_complex::Complex;

pub trait Distance<T> {
    fn dist(&self, other: &Self) -> T;
}

impl Distance<f32> for f32 {
    fn dist(&self, other: &Self) -> f32 {
        (self - other).abs()
    }
}

impl Distance<f64> for f64 {
    fn dist(&self, other: &Self) -> f64 {
        (self - other).abs()
    }
}

/// Hyperbolic distance on the upper half space
/// distance of 2 points according to Dal'bo, page 6.
#[macro_export]
macro_rules! impl_hyperbolic_distance {
    ($impl_type:ty) => {
        impl Distance<f64> for Complex<$impl_type> {
            fn dist(&self, other: &Self) -> f64 {
                if self.im <= 0.0 || other.im <= 0.0 {
                    panic!("Distance only for the hyperbolic upper half space");
                }
                let eucl_dist =
                    (self.re.dist(&other.re).powi(2) + self.im.dist(&other.im).powi(2)).sqrt();
                let x = eucl_dist / (2.0 * (self.im * other.im)).sqrt();
                // (inverse of sinh)(x) = ln(x + (x² + 1).sqrt)
                2.0 * (x + (1.0 + x.powi(2)).sqrt()).ln() as f64
            }
        }
    };
}

impl_hyperbolic_distance! {f32}
impl_hyperbolic_distance! {f64}

// TODO:  impl dist for boundary points

/// NOTE: for T in { i8, i32, i64 } etc, there is in general NO unique mid point
pub trait Mid {
    fn mid(&self, other: &Self) -> Self;
}

#[macro_export]
macro_rules! impl_mid {
    ($impl_type:ty) => {
        impl Mid for $impl_type {
            fn mid(&self, other: &Self) -> Self {
                (self + other) / 2.0
            }
        }
    };
}

impl_mid! { f32 }
impl_mid! { f64 }

/// The parametrization of a Euclidean circle in C.
pub struct EuclideanCircle<T> {
    pub center: Complex<T>,
    pub radius: T,
}

impl<T> PartialEq for EuclideanCircle<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.center == other.center && self.radius == other.radius
    }
}
impl<T> Eq for EuclideanCircle<T> where T: PartialEq {}

pub trait Drawable2d<T> {
    fn draw(&self, n_curve_points: usize) -> Vec<(T, T)>;
}

impl Drawable2d<f64> for EuclideanCircle<f64> {
    fn draw(&self, n_curve_points: usize) -> Vec<(f64, f64)> {
        let mut curve = Vec::with_capacity(n_curve_points);
        let angle_step = 2.0 * std::f64::consts::PI / (n_curve_points as f64);
        let mut angle: f64 = 0.0;
        for _ in 0..=n_curve_points {
            curve.push((
                self.center.re + self.radius * angle.cos(),
                self.center.im + self.radius * angle.sin(),
            ));
            angle += angle_step;
        }
        curve
    }
}

pub fn draw_euclidean_arc(center: f64, radius: f64, n_curve_pts: usize) -> Vec<(f64, f64)> {
    let mut curve = Vec::with_capacity(n_curve_pts);
    let angle_step = std::f64::consts::PI / (n_curve_pts as f64);
    let mut angle: f64 = 0.0;
    for _ in 0..=n_curve_pts {
        curve.push((center + radius * angle.cos(), radius * angle.sin()));
        angle += angle_step;
    }
    curve
}

// TODO: impl Mid for complex, hyperbolic points

/*/
/// The [`Busemann function`](https://en.wikipedia.org/wiki/Busemann_function) of a boundary point `$\xi$` at infinity
/// starting at the `base_point` (`$\gamma(0)$`).
pub struct BusemannParams<T> {
    pub base_point: Complex<T>,
    pub boundary_point: BoundaryPoint<T>,
}

impl<T> Default for BusemannParams<T>
where
    T: AddIdentity + MulIdentity,
{
    fn default() -> Self {
        Self {
            base_point: Complex {
                re: AddIdentity::zero(),
                im: MulIdentity::one(),
            },
            boundary_point: BoundaryPoint::Infinity,
        }
    }
}

/// The `Busemann function` of a boundary point `$\xi$` serves as a height function on the hyperbolic space:
/// `h(x) = | B(x) - B(b) |` for a fixed point `b`, or equivalently,
/// fixing a level set `$B_0$` (through `b`) via `h(x) = | B(x) - B_0 |`
/// which equals the hyperbolic distance between `x` and `B_0`.
///
/// Note that the level sets of this function correspond to horocycles based at `$\xi$` and exhaust the hyperbolic space.
pub trait Height<T> {
    fn height(&self) -> f64;
}

// should be for UpperHalfSpace
impl<T> Height<T> for Complex<T> {
    fn height(&self) -> f64 {
        1.0
    }
}

// pub trait UpperHalfSpace: SetRestriction {
//     fn condition(&self) -> bool;
// }

// impl<T> UpperHalfSpace for num_complex::Complex<T>
// where
//     T: AddIdentity,
// {
//     fn condition(&self) -> bool {
//         self.im > AddIdentity::zero()
//     }
// }



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

*/

*/
