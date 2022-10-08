use crate::{
    algebraic_extensions::{
        Group, IsPositive, MulIdentity, Numeric, NumericAddIdentity, SquareRoot,
    },
    group_action::{Action, FinitelyGeneratedGroup, SpecialLinear},
    moebius::MoebiusTransformation,
    set_extensions::{SetRestriction, Wrapper},
};
use num_complex::Complex;
use std::ops::{Add, Deref, Div, Mul, Neg};

/// Multiplicative group implementation for MoebiusTransformation satisfying determinant == 1.
impl<T> Group for MoebiusTransformation<T>
where
    Self: SpecialLinear<T> + MulIdentity + Mul<Output = MoebiusTransformation<T>> + Clone,
    T: PartialEq + Clone + Neg<Output = T>,
{
    fn combine(&self, other: &Self) -> MoebiusTransformation<T> {
        self.clone() * other.clone()
    }

    fn identity() -> Self {
        Self::one()
    }

    fn inv(&self) -> Self {
        // No need to check the determinant, since it is assumed to be 1.
        // See `MoebiusTransformation.inverse()` for the formula.
        MoebiusTransformation::<T> {
            a: self.d.clone(),
            b: -self.b.clone(),
            c: -self.c.clone(),
            d: self.a.clone(),
        }
    }
}

/// Helper:
/// We identify Moebius transformations satisfying the condition `determinant == 1` with the
/// subset [SL(2,R)](https://en.wikipedia.org/wiki/SL2(R)) within 2x2 matrices.
/// Due to mathematical limitations in rust we cannot model this condition, i.e. this subset, and hence
/// use the wrapper `SpecialLinearMoebiusTransformation<_>` which checks the condition upon construction and assumes the condition.
pub struct SpecialLinearMoebiusTransformation<T> {
    /// Moebius transformation with determinant == 1
    m: MoebiusTransformation<T>,
}

impl<T> SpecialLinearMoebiusTransformation<T> {
    fn rescale(&self, positive_determinant: T) -> Self
    where
        T: Numeric + Div<Output = T> + MulIdentity + SquareRoot + Copy,
    {
        let scalar = T::one() / positive_determinant.square_root();
        Self { m: self.m * scalar }
    }

    /// impl `TryFrom` but with a numeric threshold to check and option as result.
    pub fn try_from(m: MoebiusTransformation<T>, numeric_threshold: Option<f64>) -> Option<Self>
    where
        T: Numeric
            + NumericAddIdentity
            + Div<Output = T>
            + MulIdentity
            + SquareRoot
            + IsPositive
            + Copy,
    {
        if m.is_invertible(numeric_threshold) {
            let determinant = m.determinant();
            if determinant.is_positive() {
                let s = Self { m };
                return Some(s.rescale(determinant));
            }
        }
        None
    }
}

impl<T> Deref for SpecialLinearMoebiusTransformation<T> {
    type Target = MoebiusTransformation<T>;
    fn deref(&self) -> &Self::Target {
        &self.m
    }
}

impl<T> From<MoebiusTransformation<T>> for SpecialLinearMoebiusTransformation<T>
where
    MoebiusTransformation<T>: SpecialLinear<T>,
{
    fn from(m: MoebiusTransformation<T>) -> Self
    where
        MoebiusTransformation<T>: SpecialLinear<T>,
    {
        Self {
            m: MoebiusTransformation::try_new(m).unwrap(),
        }
    }
}

/// NOTE: `SpecialLinearMoebiusTransformation<T>` as wrapper of `MoebiusTransformation<T>` satisfying `SpecialLinear`
/// will hence implement Group under relevant trait bounds, see `algebraic_extensions.rs`
impl<T> Wrapper for SpecialLinearMoebiusTransformation<T>
where
    MoebiusTransformation<T>: SpecialLinear<T>,
{
    type Inner = MoebiusTransformation<T>;
}

impl<T> SpecialLinear<T> for SpecialLinearMoebiusTransformation<T> where
    MoebiusTransformation<T>: SpecialLinear<T>
{
}

impl<T> PartialEq for SpecialLinearMoebiusTransformation<T>
where
    MoebiusTransformation<T>: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.m == other.m
    }
}

impl<T> Eq for SpecialLinearMoebiusTransformation<T> where MoebiusTransformation<T>: PartialEq {}
impl<T> std::hash::Hash for SpecialLinearMoebiusTransformation<T>
where
    MoebiusTransformation<T>: std::hash::Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.m.hash(state);
    }
}

impl<T> Clone for SpecialLinearMoebiusTransformation<T>
where
    T: Clone,
{
    fn clone(&self) -> Self {
        Self { m: self.m.clone() }
    }
}

/// From [wikipedia](https://en.wikipedia.org/wiki/Fuchsian_group):
/// In mathematics, a Fuchsian group is a discrete subgroup of PSL(2,R).
/// The group PSL(2,R) can be regarded equivalently as a group of isometries of the hyperbolic plane, or conformal transformations of the unit disc, or conformal transformations of the upper half plane, so a Fuchsian group can be regarded as a group acting on any of these spaces.
/// There are some variations of the definition:
/// sometimes the Fuchsian group is assumed to be finitely generated,
/// sometimes it is allowed to be a subgroup of PGL(2,R) (so that it contains orientation-reversing elements), and
/// sometimes it is allowed to be a Kleinian group (a discrete subgroup of PSL(2,C)) which is conjugate to a subgroup of PSL(2,R).
///
/// NOTE: for simplicity we will in the following restrict to the case of <i>finitely generated</i> Fuchsian groups.
pub struct FuchsianGroup<T>
where
    MoebiusTransformation<T>: SpecialLinear<T>,
{
    generators: Vec<MoebiusTransformation<T>>,
    // TODO: maybe add...
    // inverse_generator: Vec<SpecialLinear<T>>,
}

impl<T> FinitelyGeneratedGroup for FuchsianGroup<T>
where
    MoebiusTransformation<T>: Group + SpecialLinear<T>,
{
    type GroupElement = MoebiusTransformation<T>;

    fn generators(&self) -> &[Self::GroupElement] {
        &self.generators
    }
}

impl<T> FuchsianGroup<T>
where
    MoebiusTransformation<T>: SpecialLinear<T>,
{
    pub fn try_push(&mut self, m: MoebiusTransformation<T>) -> bool {
        if let Some(slm) = MoebiusTransformation::try_new(m) {
            self.generators.push(slm)
        }
        false
    }

    /// Tries to create a `ProjectedMoebiusTransformation<T>` for each 'raw generator'
    /// of type `MoebiusTransformations<T>` satisfying `determinant == 1`.
    pub fn create_from_valid(raw_generators: Vec<MoebiusTransformation<T>>) -> Self {
        let mut generators = Vec::new();

        for m in raw_generators.into_iter() {
            if let Some(slm) = MoebiusTransformation::try_new(m) {
                generators.push(slm)
            }
        }

        Self { generators }
    }

    /// Tries to create a `ProjectedMoebiusTransformation<T>` for each 'raw generator',
    /// but filters for valid ones, meaning for distinct, invertible and orientation-preserving `MoebiusTransformations<T>`.
    /// For instance,
    /// - `[ -1, 0; 0, 1 ]` has determinant `-1` and is not orientation-preserving
    /// - `[ -1, 1; 0, 0 ]` has determinant `0` and is not invertible
    /// - `[ 2, 1; 1, 1 ]` and `[ 4, 2; 2, 2 ]` are projected to the same element and will result in only one generator
    pub fn create_projected(
        raw_generators: Vec<MoebiusTransformation<T>>,
        numeric_threshold: Option<f64>,
    ) -> Self
    where
        T: Numeric
            + Div<Output = T>
            + NumericAddIdentity
            + MulIdentity
            + SquareRoot
            + IsPositive
            + Copy
            + PartialEq,
        MoebiusTransformation<T>: PartialEq,
    {
        // TODO: filter out duplicates
        let generators = raw_generators
            .into_iter()
            .flat_map(|m| {
                SpecialLinearMoebiusTransformation::<T>::try_from(m, numeric_threshold)
                    .map(|slm| slm.m)
            })
            .collect::<Vec<MoebiusTransformation<T>>>();
        Self { generators }
    }
}

/// Implement Action for float types on the complex plane.
impl<T> Action<Complex<T>> for MoebiusTransformation<T>
where
    T: Numeric + Copy + PartialEq,
    Complex<T>: Div<Output = Complex<T>>,
{
    fn map(&self, x: &Complex<T>) -> Complex<T> {
        let nom = Complex::<T> {
            re: self.a * x.re + self.b,
            im: self.a * x.im,
        };
        let denom = Complex::<T> {
            re: self.c * x.re + self.d,
            im: self.c * x.im,
        };
        // TODO: check for 0?
        nom / denom
    }
}

impl<T> Action<Complex<T>> for MoebiusTransformation<Complex<T>>
where
    Complex<T>:
        Add<Output = Complex<T>> + Mul<Output = Complex<T>> + Div<Output = Complex<T>> + Clone,
{
    fn map(&self, x: &Complex<T>) -> Complex<T> {
        let nom = self.a.clone() * x.clone() + self.b.clone();
        let denom = self.c.clone() * x.clone() + self.d.clone();
        // TODO: check for 0?
        nom / denom
    }
}

#[cfg(test)]
mod tests {
    use super::FuchsianGroup;
    use crate::{
        algebraic_extensions::Group,
        fuchsian_group::SpecialLinearMoebiusTransformation,
        group_action::{Action, Orbit},
        moebius::MoebiusTransformation,
    };
    use approx::assert_abs_diff_eq;
    use num_complex::{Complex, Complex64};

    #[test]
    fn test_projection() {
        // not orentiation preserving
        let m1 = MoebiusTransformation::<f32>::new(1.0, 2.0, 3.0, 4.0);
        assert_eq!(m1.determinant(), -2.0);
        let pm1 = SpecialLinearMoebiusTransformation::try_from(m1, None);
        assert!(pm1.is_none());

        let m2 = MoebiusTransformation::<f32>::new(-1.0, -2.0, 3.0, 4.0);
        assert_eq!(m2.determinant(), 2.0);
        let pm2 = SpecialLinearMoebiusTransformation::try_from(m2, None);
        assert!(pm2.is_some());
        assert_abs_diff_eq!(pm2.unwrap().determinant(), 1.0, epsilon = f32::EPSILON);

        // numerical checks
        let m4 = MoebiusTransformation::<f32>::new(
            1.100000000002,
            2.000000000007,
            0.000000000005,
            3.000000000001,
        );
        assert_eq!(m4.determinant(), 3.3000002);
        let pm4 = SpecialLinearMoebiusTransformation::try_from(m4, None);
        assert!(pm4.is_some());
        assert_abs_diff_eq!(pm4.unwrap().determinant(), 1.0, epsilon = f32::EPSILON);

        let m5 = MoebiusTransformation::<f32>::new(0.0002, -0.0007, 0.005, 0.001);
        assert_eq!(m5.determinant(), 3.6999998e-6);
        let pm5 = SpecialLinearMoebiusTransformation::try_from(m5, None);
        assert!(pm5.is_some());
        assert_abs_diff_eq!(pm5.unwrap().determinant(), 1.0, epsilon = f32::EPSILON);
    }

    #[test]
    fn test_fuchsian_group() {
        let m1 = MoebiusTransformation::<f64>::new(1.0, 2.0, 3.0, 4.0);
        let m2 = MoebiusTransformation::<f64>::new(-1.0, -2.0, 3.0, 4.0);

        let fg = FuchsianGroup::<f64>::create_projected(vec![m1, m2], None);
        assert_eq!(fg.generators.len(), 1);
    }

    #[test]
    fn test_action_real_line() {
        let m = MoebiusTransformation::<f64>::new(1.0, 2.0, 3.0, 4.0);
        let c = Complex::new(1.0, 0.0);

        let y = m.map(&c);
        assert_eq!(y.re, 3.0 / 7.0);
        assert_eq!(y.im, 0.0);
    }

    #[test]
    fn test_action_integers() {
        let m = MoebiusTransformation::<i32>::new(0, -1, 1, 0);
        let c = Complex::new(2, 0);

        let y = m.map(&c);
        assert_eq!(y.re, 0);
        assert_eq!(y.im, 0);

        let m = MoebiusTransformation::<f32>::new(0.0, -1.0, 1.0, 0.0);
        let c = Complex::new(2.0, 0.0);

        let y = m.map(&c);
        assert_eq!(y.re, -0.5);
        assert_eq!(y.im, 0.0);
    }

    #[test]
    fn test_action_complex() {
        let m = MoebiusTransformation::<f64>::new(1.0, 2.0, 3.0, 4.0);
        let c = Complex::new(1.0, 3.0);

        let y = m.map(&c);
        assert_eq!(y.re, 48.0 / 130.0);
        assert_eq!(y.im, -6.0 / 130.0);
    }

    #[test]
    fn test_compatibility() {
        // det == 1
        let g = SpecialLinearMoebiusTransformation::try_from(
            MoebiusTransformation::<f64>::new(3.0, 2.0, 4.0, 3.0),
            None,
        );
        let h = SpecialLinearMoebiusTransformation::try_from(
            MoebiusTransformation::<f64>::new(-3.0, 2.0, -5.0, 3.0),
            None,
        );

        assert!(g.is_some());
        assert!(h.is_some());

        let g = g.unwrap();
        let h = h.unwrap();

        let mut gh = SpecialLinearMoebiusTransformation::<f64>::identity();

        let x = Complex::new(1.0, 3.0);
        let mut c = x.clone();

        for _ in 0..5 {
            gh = gh.combine(&h).combine(&g);
            c = g.map(&c);
            c = h.map(&c);
        }

        let y = gh.map(&x);

        assert_abs_diff_eq!(y.re, c.re, epsilon = 1e-15);
        assert_abs_diff_eq!(y.im, c.im, epsilon = 1e-15);
    }

    #[test]
    fn test_single_orbit() {
        // det == 1
        let m = MoebiusTransformation::<f64>::new(3.0, 2.0, 4.0, 3.0);
        let mut c = Complex::new(1.0, 3.0);

        for _ in 0..100 {
            c = m.map(&c);
            assert!(c.im >= 0.0);
        }
    }

    #[test]
    fn test_orbit() {
        let g = MoebiusTransformation::<f64>::new(3.0, 1.0, 8.0, 3.0);
        let h = MoebiusTransformation::<f64>::new(-3.0, 2.0, -5.0, 3.0);
        // horocyclic
        let j = MoebiusTransformation::<f64>::new(1.0, 10.0, 0.0, 1.0);

        let fuchsian_group = FuchsianGroup::create_from_valid(vec![g, h, j]);

        assert_eq!(fuchsian_group.generators.len(), 3);

        let base_point = Complex64::new(1.0, 1.0);
        let orbit = Orbit::create(&fuchsian_group, &base_point, 100, None);

        assert_eq!(orbit.points.len(), 100);

        for p in orbit.points {
            assert!(p.im > 0.0);
        }
    }

    #[test]
    fn test_orbit_modular_group() {
        // see https://en.wikipedia.org/wiki/Modular_group
        // the modular group is generated by the transformations z -> z+1 and z -> -1/z,
        // corresponding to the Moebius transformations g = [ 1, 1; 0, 1 ] and h = [ 0, -1; 1, 0 ]

        let g = MoebiusTransformation::<i32>::new(1, 1, 0, 1);
        let h = MoebiusTransformation::<i32>::new(0, -1, 1, 0);

        let fuchsian_group = FuchsianGroup::create_from_valid(vec![g, h]);

        assert_eq!(fuchsian_group.generators.len(), 2);

        let base_point = Complex::<i32>::new(1, 0);
        let orbit = Orbit::create(&fuchsian_group, &base_point, 100, None);

        assert_eq!(orbit.points.len(), 100);

        for p in orbit.points {
            assert_eq!(p.im, 0);
        }

        let g = MoebiusTransformation::<f32>::new(1.0, 1.0, 0.0, 1.0);
        let h = MoebiusTransformation::<f32>::new(0.0, -1.0, 1.0, 0.0);

        let fuchsian_group = FuchsianGroup::create_from_valid(vec![g, h]);

        assert_eq!(fuchsian_group.generators.len(), 2);

        let base_point = Complex::<f32>::new(1.0, 0.0);
        let orbit = Orbit::create(&fuchsian_group, &base_point, 100, None);

        assert_eq!(orbit.points.len(), 100);

        for p in orbit.points {
            assert_eq!(p.im, 0.0);
            assert_ne!(p.re, 0.0);
        }
    }
}
