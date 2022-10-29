#![warn(dead_code)]

pub(crate) mod algebraic_extensions;
pub mod fuchsian_group;
pub mod geometry;
pub mod group_action;
pub mod moebius;
pub(crate) mod set_extensions;

pub const NUMERIC_THRESHOLD: f64 = 1e-16;
