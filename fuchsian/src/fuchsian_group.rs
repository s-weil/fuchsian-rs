use crate::{
    algebraic_extensions::{MulIdentityElement, Numeric, NumericAddIdentity},
    moebius::MoebiusTransformation,
};
use std::{
    collections::HashSet,
    ops::{Deref, Div},
};

/// The [mathematical group](https://en.wikipedia.org/wiki/Group_(mathematics)#Definition)
/// definition except for the associativity identity.
/// In particular, `MoebiusTransformation<T>` is an additive group
/// ```
/// impl<T> Group for T
/// where
///     T: AddIdentityElement + Add<Output = Self> + Neg<Output = Self> + Copy + Eq,
/// {
///     fn combine(&self, other: &Self) -> Self {
///         *self + *other
///     }
///
///     fn identity() -> Self {
///         Self::zero()
///     }
///
///     fn inverse(&self) -> Self {
///         -*self
///     }
/// }
/// ```
pub trait Group {
    /// The binary operation.
    fn combine(&self, other: &Self) -> Self;

    /// The identity element.
    fn identity() -> Self;

    /// The inverse element.
    fn inverse(&self) -> Self;

    // how to model associativity identity ?
    // fn swap((a, b) : (Self, Self), c: Self) -> (Self, (Self, Self));
    // fn is_associative(&self, b: &Self, c: &Self) -> bool {
    //     let ab_c = self.combine(b).combine(c);
    //     let a_bc = self.combine(b.combine(c));
    //     ab_c == ab_c
    // }
}

/// Helper
/// We would like to model MoebiusTransformation with the restriction of determinant == 1.
/// Not possible due to mathematical limitations in rust.
// type ProjectedMoebiusTransformation<T> = MoebiusTransformation<T>;

// enum ProjectedMoebiusTransformation2<'a, T> {
//     /// A MoebiusTransformation with the restriction of determinant == 1
//     SpecialLinear(&'a MoebiusTransformation<T>),
//     Projected(MoebiusTransformation<T>),
// }

/// Helper:
/// We would like to model MoebiusTransformation with the condition `determinant == 1`.
/// Not possible due to mathematical limitations in rust.
/// So we have the wrapper `ProjectedMoebiusTransformation<_>` which checks the condition upon construction.
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
        // *self * *other
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

/*
From [wikipedia](https://en.wikipedia.org/wiki/Fuchsian_group):
In mathematics, a Fuchsian group is a discrete subgroup of PSL(2,R).
The group PSL(2,R) can be regarded equivalently as a group of isometries of the hyperbolic plane, or conformal transformations of the unit disc, or conformal transformations of the upper half plane, so a Fuchsian group can be regarded as a group acting on any of these spaces.
There are some variations of the definition:
sometimes the Fuchsian group is assumed to be finitely generated,
sometimes it is allowed to be a subgroup of PGL(2,R) (so that it contains orientation-reversing elements), and
sometimes it is allowed to be a Kleinian group (a discrete subgroup of PSL(2,C)) which is conjugate to a subgroup of PSL(2,R).

NOTE: in the following we will restrict to the case of <i>finitely generated</i> Fuchsian groups
*/
pub trait FinitelyGeneratedGroup: Group + Sized {
    fn generators() -> HashSet<Self>;
}

pub struct ProjectedMoebiusTransformationGroup<T> {
    generators: HashSet<ProjectedMoebiusTransformation<T>>,
}

impl<T> ProjectedMoebiusTransformationGroup<T> {
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

// struct ProjectedMoebiusTransformation2<T> {
//     m: MoebiusTransformation<T>,
// }

// impl<T> Deref for ProjectedMoebiusTransformation2<T> {
//     type Target = MoebiusTransformation<T>;

//     fn deref(&self) -> &Self::Target {
//         &self.m
//     }
// }

// impl<T> ProjectedMoebiusTransformation2<T> {
//     fn try_project(m: MoebiusTransformation<T>, threshold: Option<f64>) -> Option<Self>
//     where
//         T: Numeric + NumericAddIdentity + Div<Output = T> + MulIdentityElement + std::marker::Copy,
//     {
//         if m.is_invertible(threshold) {
//             let determinant = m.determinant();

//             let scalar = T::one() / (determinant * determinant); // TODO: need signed
//             return Some(ProjectedMoebiusTransformation2 { m: m * scalar });
//         }
//         None
//     }
// }

// /// Group implementation for additive groups
// /// In particular, MoebiusTransformation<_> is an additive group
// // impl<T> Group for T
// // where
// //     T: AddIdentityElement + Add<Output = Self> + Neg<Output = Self> + Copy + Eq,
// // {
// //     fn combine(&self, other: &Self) -> Self {
// //         *self + *other
// //     }

// //     fn identity() -> Self {
// //         Self::zero()
// //     }

// //     fn inverse(&self) -> Self {
// //         -*self
// //     }
// // }

// pub struct discrete Group
// trait GroupAction on disk / upper half space
// orbit
//
