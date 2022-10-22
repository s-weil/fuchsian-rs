use fuchsian::geometry::boundary::BoundaryPoint;
use fuchsian::geometry::geodesics::{Arc, GeodesicBoundary, GeodesicLine};
use fuchsian::geometry::horocycle::{EuclieanCircle, HoroCycle, ParametricHoroCycle};
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
    // TODO: check for SL
    MoebiusTransformation::new(a, b, c, d)
}

fn parse_complex((re, im): (f64, f64)) -> Complex<f64> {
    Complex::new(re, im)
}

fn parse_geodesic_boundary((start, end): (f64, f64)) -> GeodesicBoundary<f64> {
    GeodesicBoundary::new(BoundaryPoint::Regular(start), BoundaryPoint::Regular(end))
}

fn parse_horocycle(height: f64) -> HoroCycle<f64> {
    HoroCycle::new(BoundaryPoint::Infinity, height)
}

fn plot_euclidean_arc(arc: Arc<f64>, n_curve_pts: usize) -> Vec<(f64, f64)> {
    let mut curve = Vec::with_capacity(n_curve_pts);
    let angle_step = std::f64::consts::PI / (n_curve_pts as f64);
    let mut angle: f64 = 0.0;
    for _ in 0..=n_curve_pts {
        curve.push((
            arc.center + arc.radius * angle.cos(),
            arc.radius * angle.sin(),
        ));
        angle += angle_step;
    }
    curve
}

fn plot_geodesic(geodesic_boundary: GeodesicBoundary<f64>, n_curve_pts: usize) -> Vec<(f64, f64)> {
    let line = GeodesicLine::from(geodesic_boundary);

    match line {
        GeodesicLine::Line(b) => {
            vec![(b, 0.0), (b, n_curve_pts as f64)]
        }
        GeodesicLine::Arc(arc) => plot_euclidean_arc(arc, n_curve_pts),
    }
}

fn plot_euclidean_circle(circle: EuclieanCircle<f64>, n_curve_pts: usize) -> Vec<(f64, f64)> {
    let mut curve = Vec::with_capacity(n_curve_pts);
    let angle_step = 2.0 * std::f64::consts::PI / (n_curve_pts as f64);
    let mut angle: f64 = 0.0;
    for _ in 0..=n_curve_pts {
        curve.push((
            circle.center.re + circle.radius * angle.cos(),
            circle.center.im + circle.radius * angle.sin(),
        ));
        angle += angle_step;
    }
    curve
}

fn plot_horocycle(horocycle: HoroCycle<f64>, n_curve_pts: usize) -> Vec<(f64, f64)> {
    let eucl_hc = ParametricHoroCycle::from(horocycle);

    match eucl_hc {
        ParametricHoroCycle::Line(height) => {
            vec![
                (-(n_curve_pts as f64), height),
                (n_curve_pts as f64, height),
            ]
        }
        ParametricHoroCycle::Circle(eucl_circle) => plot_euclidean_circle(eucl_circle, n_curve_pts),
    }
}

fn parse_generator_mode(mode: Option<String>) -> Option<PickGeneratorMode> {
    if let Some(m) = mode {
        match m.to_lowercase().trim() {
            "random" => Some(PickGeneratorMode::Random),
            "sequential" => Some(PickGeneratorMode::Sequential),
            _ => todo!(), // return Err(PyErr::new("Invalid mode")),
        }
    } else {
        None
    }
}

#[pyfunction]
fn geodesic_orbit(
    moebius: Vec<((f64, f64), (f64, f64))>,
    end_points: (f64, f64),
    n_pts: usize,
    n_curve_pts: usize,
    mode: Option<String>,
) -> PyResult<Vec<Vec<(f64, f64)>>> {
    let m = moebius.into_iter().map(parse_moebius).collect();
    let base_geodesic = parse_geodesic_boundary(end_points);
    let pick_mode = parse_generator_mode(mode);

    let fuchsian_group = FuchsianGroup::create_projected(m, None);
    let orbit = Orbit::sample(&fuchsian_group, &base_geodesic, n_pts, pick_mode);

    let orbit_ves = orbit
        .points
        .into_iter()
        .map(|gb| plot_geodesic(gb, n_curve_pts))
        .collect();

    Ok(orbit_ves)
}

#[pyfunction]
fn horocyclic_orbit(
    moebius: Vec<((f64, f64), (f64, f64))>,
    height: f64,
    n_pts: usize,
    n_curve_pts: usize,
    mode: Option<String>,
) -> PyResult<Vec<Vec<(f64, f64)>>> {
    let m = moebius.into_iter().map(parse_moebius).collect();
    let base_horocycle = parse_horocycle(height);
    let pick_mode = parse_generator_mode(mode);

    let fuchsian_group = FuchsianGroup::create_projected(m, None);
    let orbit = Orbit::sample(&fuchsian_group, &base_horocycle, n_pts, pick_mode);

    let orbit_ves = orbit
        .points
        .into_iter()
        .map(|gb| plot_horocycle(gb, n_curve_pts))
        .collect();

    Ok(orbit_ves)
}

#[pyfunction]
fn orbit(
    moebius: Vec<((f64, f64), (f64, f64))>,
    base_point: (f64, f64),
    n_pts: usize,
    mode: Option<String>,
) -> PyResult<Vec<(f64, f64)>> {
    let m = moebius.into_iter().map(parse_moebius).collect();
    let base_point = parse_complex(base_point);
    let pick_mode = parse_generator_mode(mode);

    let fuchsian_group = FuchsianGroup::create_projected(m, None);
    let orbit = Orbit::sample(&fuchsian_group, &base_point, n_pts, pick_mode);

    let orbit_ves = orbit.points.into_iter().map(|c| (c.re, c.im)).collect();

    Ok(orbit_ves)
}

/// A Python module implemented in Rust.
#[pymodule]
fn python_api(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(moebius_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(orbit, m)?)?;
    m.add_function(wrap_pyfunction!(geodesic_orbit, m)?)?;
    m.add_function(wrap_pyfunction!(horocyclic_orbit, m)?)?;
    Ok(())
}