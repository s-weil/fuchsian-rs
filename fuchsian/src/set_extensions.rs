/// Model a restriction to a set.
/// For a set $S$, we want to model a restriction $r$ to a subset
/// ${ x \in S : r(x) }$
pub trait Subset {
    fn restriction(&self) -> bool;

    // TODO: better modelling?
    // TODO: maybe return Error in case restriction check fails
    fn check_restriction(&self) {
        assert!(self.restriction());
    }
}
