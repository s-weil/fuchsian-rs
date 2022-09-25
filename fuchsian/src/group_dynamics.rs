use std::collections::HashSet;

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

pub trait FinitelyGeneratedGroup<M>
where
    M: Group + std::hash::Hash,
{
    fn generators(&self) -> &HashSet<M>;
}

/// The [mathematical group action](https://en.wikipedia.org/wiki/Group_action)
pub trait GroupAction {
    /// The group operating on the space
    type OpGroup: Group;
    type Space;

    fn action(g: &Self::OpGroup, x: &Self::Space) -> Self::Space;
}

// TODO:
// - Orbit
// - Fundamental Domain
