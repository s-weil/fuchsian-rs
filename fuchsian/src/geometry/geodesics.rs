pub enum Geodesic<T> {
    /// A half-circle with center on the real axis within C
    Arc(Arc<T>),
    /// A half-line perpendicular to the real axis within C
    Line(T), // touchpoint
}

pub struct Arc<T> {
    center: T,
    radius: T,
}
