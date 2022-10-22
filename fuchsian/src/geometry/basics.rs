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

// TODO: use macro instead! impl dist for complex, hyperbolic, boundary points

/// NOTE: for T in { i8, i32, i64 } etc, there is in general NO unique mid point
pub trait Mid {
    fn mid(&self, other: &Self) -> Self;
}

impl Mid for f32 {
    fn mid(&self, other: &Self) -> Self {
        (self + other) / 2.0
    }
}

impl Mid for f64 {
    fn mid(&self, other: &Self) -> Self {
        (self + other) / 2.0
    }
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
*/
