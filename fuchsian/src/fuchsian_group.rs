use crate::{
    algebraic_extensions::{MulIdentityElement, Numeric, NumericAddIdentity},
    group_dynamics::{FinitelyGeneratedGroup, Group},
    moebius::MoebiusTransformation,
};
use std::{
    collections::HashSet,
    ops::{Deref, Div},
};

/// Helper:
/// We would like to model MoebiusTransformation with the condition `determinant == 1`,
/// which is not possible due to mathematical limitations in rust.
/// Use the wrapper `ProjectedMoebiusTransformation<_>` which checks the condition upon construction.
struct ProjectedMoebiusTransformation<T> {
    m: MoebiusTransformation<T>,
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

impl<T> ProjectedMoebiusTransformation<T> {
    fn rescale(&self, determinant: T) -> Self
    where
        T: Numeric + Div<Output = T> + MulIdentityElement + Copy,
    {
        let scalar = T::one() / (determinant * determinant); // TODO: need to be signed
        Self { m: self.m * scalar }
    }

    /// impl `TryFrom` but with a numeric threshold to check.
    pub fn try_from(m: MoebiusTransformation<T>, numeric_threshold: Option<f64>) -> Option<Self>
    where
        T: Numeric + NumericAddIdentity + Div<Output = T> + MulIdentityElement + Copy,
    {
        if m.is_invertible(numeric_threshold) {
            let determinant = m.determinant();
            let s = Self { m };
            return Some(s.rescale(determinant));
        }
        None
    }
}

// Multiplicative group implementation for MoebiusTransformation with the restriction of determinant == 1.
impl<T> Group for ProjectedMoebiusTransformation<T>
where
    T: Numeric + MulIdentityElement + Copy,
{
    fn combine(&self, _other: &Self) -> Self {
        let m1 = self.m;
        let m2 = self.m;
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
    generators: HashSet<ProjectedMoebiusTransformation<T>>,
}

impl<T> FinitelyGeneratedGroup<ProjectedMoebiusTransformation<T>> for FuchsianGroup<T>
where
    T: Numeric + MulIdentityElement + Copy,
    MoebiusTransformation<T>: PartialEq + std::hash::Hash,
{
    fn generators(&self) -> &HashSet<ProjectedMoebiusTransformation<T>> {
        &self.generators
    }
}

impl<T> FuchsianGroup<T> {
    pub fn for_valid(raw_generators: Vec<MoebiusTransformation<T>>) -> Self
    where
        T: Numeric + MulIdentityElement + Eq + Copy,
        MoebiusTransformation<T>: PartialEq + std::hash::Hash,
    {
        let mut generators = HashSet::new();

        let one = T::one();
        for m in raw_generators.into_iter() {
            if m.determinant() == one {
                generators.insert(ProjectedMoebiusTransformation { m });
            }
        }

        Self { generators }
    }

    pub fn create_projected(
        raw_generators: Vec<MoebiusTransformation<T>>,
        numeric_threshold: Option<f64>,
    ) -> Self
    where
        T: Numeric
            + Div<Output = T>
            + NumericAddIdentity
            + MulIdentityElement
            + std::marker::Copy
            + PartialEq,
        MoebiusTransformation<T>: PartialEq + std::hash::Hash,
    {
        let generators = raw_generators
            .into_iter()
            .flat_map(|m| ProjectedMoebiusTransformation::<T>::try_from(m, numeric_threshold))
            .collect::<HashSet<ProjectedMoebiusTransformation<T>>>();
        Self { generators }
    }
}
