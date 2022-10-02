use crate::{
    algebraic_extensions::{IsPositive, MulIdentity, Numeric, NumericAddIdentity, SquareRoot},
    group_dynamics::{Action, FinitelyGeneratedGroup, Group},
    moebius::MoebiusTransformation,
};
use num_complex::Complex;
use std::ops::{Deref, Div};

/// Helper:
/// We identify MoebiusTransformations with the condition `determinant == 1` with the orientation preserving subset PSL(2,R) within 2x2 matrices.
/// Due to mathematical limitations in rust we cannot model this condition, i.e. this subset, and hence
/// use the wrapper `ProjectedMoebiusTransformation<_>` which checks the condition upon construction and assumes the condition.
struct ProjectedMoebiusTransformation<T> {
    pub(crate) m: MoebiusTransformation<T>,
}

impl<T> Deref for ProjectedMoebiusTransformation<T> {
    type Target = MoebiusTransformation<T>;
    fn deref(&self) -> &Self::Target {
        &self.m
    }
}

impl<T> PartialEq for ProjectedMoebiusTransformation<T>
where
    MoebiusTransformation<T>: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.m == other.m
    }
}

impl<T> Eq for ProjectedMoebiusTransformation<T> where MoebiusTransformation<T>: PartialEq {}
impl<T> std::hash::Hash for ProjectedMoebiusTransformation<T>
where
    MoebiusTransformation<T>: std::hash::Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.m.hash(state);
    }
}

// TODO: revisit the trait bounds below; use only what is actually required!

impl<T> ProjectedMoebiusTransformation<T> {
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

// Multiplicative group implementation for MoebiusTransformation with the restriction of determinant == 1.
impl<T> Group for ProjectedMoebiusTransformation<T>
where
    T: Numeric + MulIdentity + Copy + Eq,
{
    fn combine(&self, other: &Self) -> Self {
        let m1 = self.m;
        let m2 = other.m;
        Self { m: m1 * m2 }
    }

    fn identity() -> Self {
        let one = MoebiusTransformation::one();
        Self { m: one }
    }

    fn inverse(&self) -> Self {
        // Do not check the determinant, it is assumed to be 1.
        // Check `MoebiusTransformation.inverse()` for the formula.
        let inverse = MoebiusTransformation::<T> {
            a: *&self.m.d,
            b: -*&self.m.b,
            c: -*&self.m.c,
            d: *&self.m.a,
        };
        Self { m: inverse }
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
pub struct FuchsianGroup<T> {
    generators: Vec<ProjectedMoebiusTransformation<T>>,
}

impl<T> FinitelyGeneratedGroup<ProjectedMoebiusTransformation<T>> for FuchsianGroup<T>
where
    ProjectedMoebiusTransformation<T>: Group, // T: Numeric + MulIdentityElement + Copy + Eq,
                                              // MoebiusTransformation<T>: PartialEq + std::hash::Hash,
{
    fn generators(&self) -> &[ProjectedMoebiusTransformation<T>] {
        &self.generators
    }
}

impl<T> FuchsianGroup<T> {
    /// Tries to create a `ProjectedMoebiusTransformation<T>` for each 'raw generator',
    /// but filters for distinct `MoebiusTransformations<T>` with `determinant == 1`.
    pub fn create_valid(raw_generators: Vec<MoebiusTransformation<T>>) -> Self
    where
        T: Numeric + MulIdentity + Eq + Copy,
        MoebiusTransformation<T>: PartialEq,
    {
        let mut generators = Vec::new();

        let one = T::one();
        for m in raw_generators.into_iter() {
            if m.determinant() == one {
                generators.push(ProjectedMoebiusTransformation { m });
            }
        }

        Self { generators }
    }

    /// Tries to create a `ProjectedMoebiusTransformation<T>` for each 'raw generator',
    /// but filters for valid ones, meaning for distinct, invertible and orientation-preserving `MoebiusTransformations<T>`.
    /// For instance,
    /// - `[ -1, 0; 0, 1 ]` has determinant `-1` and is not orientation-preserving
    /// - `[ -1, 1; 0, 0 ]` has determinant `0` and is not invertible
    /// - `[ 2, 1; 1, 1 ]` and `[ 4, 2; 2, 2 ]` are projected to the same element and will result in only one generator // TODO:
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
        let generators = raw_generators
            .into_iter()
            .flat_map(|m| ProjectedMoebiusTransformation::<T>::try_from(m, numeric_threshold))
            .collect::<Vec<ProjectedMoebiusTransformation<T>>>();
        Self { generators }
    }
}

/// Implement Action for float types on the complex plane.
impl<T> Action for ProjectedMoebiusTransformation<T>
where
    T: Eq + Numeric + Copy,
    Complex<T>: Div<Output = Complex<T>>,
{
    type Space = num_complex::Complex<T>;

    fn action(&self, x: &Self::Space) -> Self::Space {
        let nom = Complex::<T> {
            re: self.m.a * x.re + self.m.b,
            im: self.m.a * x.im,
        };
        let denom = Complex::<T> {
            re: self.m.c * x.re + self.m.d,
            im: self.m.c * x.im,
        };
        // TODO: check for 0?
        nom / denom
    }
}

// impl<T> GroupAction for FuchsianGroup<T>
// where
//     T: Eq + Sized,
// {
//     type Space = num_complex::Complex<T>;

//     fn action(&self,  x: &Self::Space) -> Self::Space {
//         self.
//     }

//     // how to model identities in general?
//     fn identity_check(&self, x: Self::Space) -> bool {
//         self.action(&Group::identity(), &x) == x
//     }

//     fn compatibility_check(&self, g: &Self::OpGroup, h: &Self::OpGroup, x: &Self::Space) -> bool {
//         let g_hx = self.action(g, &self.action(h, x));
//         let gh_x = self.action(&g.combine(h), x);
//         g_hx == gh_x
//     }
// }

#[cfg(test)]
mod tests {
    use crate::{fuchsian_group::ProjectedMoebiusTransformation, moebius::MoebiusTransformation};
    use approx::assert_abs_diff_eq;

    use super::FuchsianGroup;

    #[test]
    fn test_projection() {
        // not orentiation preserving
        let m1 = MoebiusTransformation::<f32>::new(1.0, 2.0, 3.0, 4.0);
        assert_eq!(m1.determinant(), -2.0);
        let pm1 = ProjectedMoebiusTransformation::try_from(m1, None);
        assert!(pm1.is_none());

        let m2 = MoebiusTransformation::<f32>::new(-1.0, -2.0, 3.0, 4.0);
        assert_eq!(m2.determinant(), 2.0);
        let pm2 = ProjectedMoebiusTransformation::try_from(m2, None);
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
        let pm4 = ProjectedMoebiusTransformation::try_from(m4, None);
        assert!(pm4.is_some());
        assert_abs_diff_eq!(pm4.unwrap().determinant(), 1.0, epsilon = f32::EPSILON);

        let m5 = MoebiusTransformation::<f32>::new(0.0002, -0.0007, 0.005, 0.001);
        assert_eq!(m5.determinant(), 3.6999998e-6);
        let pm5 = ProjectedMoebiusTransformation::try_from(m5, None);
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
}
