/// Model a restriction to a set.
/// For a set $S$, we want to model a restriction (or property) $r$ to a subset
/// $S_sub = { x \in S : r(x) = true }$
pub trait SetRestriction: Sized {
    fn restriction(s: &Self) -> bool;

    // TODO: consider Error instead of option
    /// This constructor for $S_sub$ shall be used to make sure the restriction holds.
    fn try_new(s: Self) -> Option<Self> {
        if Self::restriction(&s) {
            Some(s)
        } else {
            None
        }
    }
}
