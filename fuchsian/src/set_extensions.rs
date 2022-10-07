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

pub trait Wrapper<S>: From<S> + Deref<Target = S> {}
