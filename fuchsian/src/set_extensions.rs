use std::ops::Deref;

/// Model a restriction to a set.
/// For a set $S$, we want to model a restriction (i.e. a holding condition) $c$ to a subset
/// $S_sub = { x \in S : c(x) = true }$
/// Due to the mathematical limitations, we require each `constructed` element to satisfy the condition.
pub trait SetRestriction: Sized {
    fn condition(&self) -> bool;

    // TODO: consider Error instead of option
    /// This constructor for $S_sub$ makes sure the condition holds if created.
    fn try_new(s: Self) -> Option<Self> {
        if s.condition() {
            Some(s)
        } else {
            None
        }
    }
}

/// Wrapper holding an inner type `Inner`.
/// NOTE: Although the `Deref<Target = Inner>` is sufficient for many auto trait implementations,
///       we will always require full `Wrapper` to make sure it can be converted both ways.
pub trait Wrapper: From<Self::Inner> + Deref<Target = Self::Inner> {
    type Inner;
}

impl<M, S> SetRestriction for M
where
    M: Wrapper<Inner = S>,
    S: SetRestriction,
{
    fn condition(&self) -> bool {
        self.deref().condition()
    }
}
