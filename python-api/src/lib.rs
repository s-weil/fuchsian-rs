use pyo3::prelude::*;
extern crate fuchsian;
use fuchsian::moebius::MoebiusTransformation;

// TODO: rename project into pyperbolic, or something with python and fuchsian / hyperbolic

// run `maturin develop` after changes

#[pyfunction]
fn moebius(a: f32, b: f32, c: f32, d: f32) -> PyResult<Vec<f32>> {
    let m = MoebiusTransformation::new(a, b, c, d);
    let mut v = Vec::with_capacity(4);
    v.push(m.a);
    v.push(m.b);
    v.push(m.c);
    v.push(m.d);

    Ok(v)
    // Python::with_gil(|py| m.to_object(py));
}

/// A Python module implemented in Rust.
#[pymodule]
fn python_api(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(moebius, m)?)?;
    Ok(())
}
