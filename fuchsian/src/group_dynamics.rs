// use crate::algebraic_extensions::AddIdentity;
// // The [mathematical group](https://en.wikipedia.org/wiki/Group_(mathematics)#Definition)
// // definition except for the associativity identity.
// // /// In particular, `MoebiusTransformation<T>` is an additive group
// /// ```
// /// use std::ops::{Add, Neg};
// /// impl<T> Group for T
// /// where
// ///     T: AddIdentityElement + Add<Output = Self> + Neg<Output = Self> + Copy + PartialEq + Sized,
// /// {
// ///     fn combine(&self, other: &Self) -> Self {
// ///         *self + *other
// ///     }
// ///
// ///     fn identity() -> Self {
// ///         Self::zero()
// ///     }
// ///
// ///     fn inverse(&self) -> Self {
// ///         -*self
// ///     }
// /// }
// /// ```
// pub trait Group: PartialEq + Sized {
//     /// The binary operation.
//     fn combine(&self, other: &Self) -> Self;

//     /// The identity element.
//     fn identity() -> Self;

//     /// The inverse element.
//     fn inverse(&self) -> Self;

//     // how to model identities in general?
//     fn associativity_check(&self, b: &Self, c: &Self) -> bool {
//         let ab_c = self.combine(b).combine(c);
//         let a_bc = self.combine(&b.combine(c));
//         ab_c == a_bc
//     }
// }

use crate::algebraic_extensions::Group;

pub trait FinitelyGeneratedGroup<M>
where
    M: Group,
{
    fn generators(&self) -> &[M];
}

pub trait Action<Space> {
    fn map(&self, x: &Space) -> Space;
}

/// The mathematical (left) [group action](https://en.wikipedia.org/wiki/Group_action) of
/// a group `G` acting on a space `X`.
pub trait GroupAction<Space>: Action<Space> + Group
where
    Space: PartialEq,
{
    // how to model identities in general?
    fn identity_check(&self, x: &Space) -> bool {
        ((&Group::identity() as &Self).map(&x)).eq(&x)
    }

    fn compatibility_check(&self, h: &Self, x: &Space) -> bool {
        let g_hx = self.map(&h.map(&x));
        let gh_x = (self.combine(h)).map(x); //.action(&g.combine(h), x);
        g_hx.eq(&gh_x)
    }
}

// TODO:
// - revisit trait definitions and required bounds
// - Orbit
// - Fundamental Domain
// - check performance of group action / orbit. i.e. what is faster, combining group elements and then action, or iterative actions
