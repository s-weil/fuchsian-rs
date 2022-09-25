use std::collections::HashSet;

/// The [mathematical group](https://en.wikipedia.org/wiki/Group_(mathematics)#Definition)
/// definition except for the associativity identity.
// /// In particular, `MoebiusTransformation<T>` is an additive group
// /// ```
// /// use std::ops::{Add, Neg};
// /// use crate::algebraic_extensions::AddIdentityElement;
// /// impl<T> Group for T
// /// where
// ///     T: AddIdentityElement + Add<Output = Self> + Neg<Output = Self> + Copy + Eq + Sized,
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
pub trait Group: Eq + Sized {
    /// The binary operation.
    fn combine(&self, other: &Self) -> Self;

    /// The identity element.
    fn identity() -> Self;

    /// The inverse element.
    fn inverse(&self) -> Self;

    // how to model identities in general?
    fn associativity_check(&self, b: &Self, c: &Self) -> bool {
        let ab_c = self.combine(b).combine(c);
        let a_bc = self.combine(&b.combine(c));
        ab_c == a_bc
    }
}

pub trait FinitelyGeneratedGroup<M>
where
    M: Group + std::hash::Hash,
{
    fn generators(&self) -> &HashSet<M>;
}

/// The mathematical (left) [group action](https://en.wikipedia.org/wiki/Group_action) of
/// a group `G` acting on a space `X`.
pub trait GroupAction {
    /// The group operating on the space
    type OpGroup: Group;
    type Space: Eq;

    fn action(&self, g: &Self::OpGroup, x: &Self::Space) -> Self::Space;

    // how to model identities in general?
    fn identity_check(&self, x: Self::Space) -> bool {
        self.action(&Group::identity(), &x) == x
    }

    fn compatibility_check(&self, g: &Self::OpGroup, h: &Self::OpGroup, x: &Self::Space) -> bool {
        let g_hx = self.action(g, &self.action(h, x));
        let gh_x = self.action(&g.combine(h), x);
        g_hx == gh_x
    }
}

// TODO:
// - revisit trait definitions and required bounds
// - Orbit
// - Fundamental Domain
// - check performance of group action / orbit. i.e. what is faster, combining group elements and then action, or iterative actions
