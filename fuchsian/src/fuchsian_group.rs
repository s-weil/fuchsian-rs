use crate::{
    algebraic_extensions::{AddIdentityElement, MulIdentityElement, Numeric, NumericAddIdentity},
    moebius::MoebiusTransformation,
};
use std::{
    collections::HashSet,
    ops::{Add, Div, Neg},
};

// pub struct discrete Group
// trait GroupAction on disk / upper half space
// orbit
//

/// See https://en.wikipedia.org/wiki/Group_(mathematics)#Definition
pub trait Group {
    /// The binary operation.
    fn combine(&self, other: &Self) -> Self;

    /// The identity element.
    fn identity() -> Self;

    /// The inverse element.
    fn inverse(&self) -> Self;

    // how to model associativity ?
    // fn swap((a, b) : (Self, Self), c: Self) -> (Self, (Self, Self));
    // fn check_associativity(&self, b: &Self, c: &Self) -> bool {
    //     let ab_c = self.combine(b).combine(c);
    //     let a_bc = self.combine(b.combine(c));
    //     ab_c.combine(&a_bc.inverse()) // eq Self::identity
    // }
}

/// Group implementation for additive groups
/// In particular, MoebiusTransformation<_> is an additive group
impl<T> Group for T
where
    T: AddIdentityElement + Add<Output = Self> + Neg<Output = Self> + Copy,
{
    fn combine(&self, other: &Self) -> Self {
        *self + *other
    }

    fn identity() -> Self {
        Self::zero()
    }

    fn inverse(&self) -> Self {
        -*self
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

pub struct ProjectedMoebiusTransformation<T> {
    pub generators: HashSet<MoebiusTransformation<T>>,
}

impl<T> ProjectedMoebiusTransformation<T> {
    fn try_project(
        m: MoebiusTransformation<T>,
        threshold: Option<f64>,
    ) -> Option<MoebiusTransformation<T>>
    where
        T: Numeric + NumericAddIdentity + Div<Output = T> + MulIdentityElement + std::marker::Copy,
    {
        if m.is_invertible(threshold) {
            let determinant = m.determinant();

            let scalar = T::one() / (determinant * determinant); // TODO: need signed
            return Some(m * scalar);
        }
        None
    }

    pub fn new(raw_generators: Vec<MoebiusTransformation<T>>, threshold: Option<f64>) -> Self
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
            .flat_map(|m| Self::try_project(m, threshold))
            .collect::<HashSet<MoebiusTransformation<T>>>();
        Self { generators }
    }
}
