use pyo3::prelude::*;
extern crate fuchsian;
use fuchsian::group_action::{Orbit, PickGeneratorMode};
use fuchsian::{fuchsian_group::FuchsianGroup, moebius::MoebiusTransformation};
use num_complex::Complex;
// TODO: rename project into pyperbolic, or something with python and fuchsian / hyperbolic

// run `maturin develop` after changes
// read https://pyo3.rs/v0.2.7/overview.html

#[pyfunction]
fn moebius_matrix(a: f64, b: f64, c: f64, d: f64) -> PyResult<((f64, f64), (f64, f64))> {
    let m = MoebiusTransformation::new(a, b, c, d);

    Ok(((m.a, m.b), (m.c, m.d)))
    // Python::with_gil(|py| m.to_object(py));
}

fn parse_moebius(((a, b), (c, d)): ((f64, f64), (f64, f64))) -> MoebiusTransformation<f64> {
    // TODO: checks

    let m = MoebiusTransformation::new(a, b, c, d);
    m
}

fn parse_complex((re, im): (f64, f64)) -> Complex<f64> {
    // TODO: checks
    let c = Complex::new(re, im);
    c
}

#[pyfunction]
fn orbit(
    moebius: Vec<((f64, f64), (f64, f64))>,
    base_point: (f64, f64),
    n_pts: usize,
    mode: Option<String>,
) -> PyResult<Vec<(f64, f64)>> {
    let m = moebius.into_iter().map(|m| parse_moebius(m)).collect();
    let base_point = parse_complex(base_point);
    let pick_mode = if let Some(m) = mode {
        match m.to_lowercase().trim() {
            "random" => Some(PickGeneratorMode::Random),
            "sequential" => Some(PickGeneratorMode::Sequential),
            _ => todo!(), // return Err(PyErr::new("Invalid mode")),
        }
    } else {
        None
    };

    let fuchsian_group = FuchsianGroup::create_projected(m, None);
    let orbit = Orbit::create(&fuchsian_group, &base_point, n_pts, pick_mode);

    let orbit_ves = orbit.points.into_iter().map(|c| (c.re, c.im)).collect();

    Ok(orbit_ves)
}

/// A Python module implemented in Rust.
#[pymodule]
fn python_api(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(moebius_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(orbit, m)?)?;
    Ok(())
}
