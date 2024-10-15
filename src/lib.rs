use std::f64::consts::PI;

use wolfram_library_link::{export, NumericArray};

mod billiard;
use billiard::Billiard;

#[export]
fn minimize(si: f64, sf: f64, n: i64, delta: f64) -> NumericArray<f64> {
    //let z = |s: &f64| [*s, (PI * *s).sin()];
    //let z = |s: &f64|  [(1.25 + s.cos()) * s.cos(), (1.25 + s.cos()) * s.sin()];
    let z = |s: &f64|  [2.0 * s.cos(), 0.5 * s.sin()];
    //let z = |s: &f64|  [s.cos(), s.sin()];
    let bill = Billiard::new(si, sf, n, delta, z);
    vec_to_na(bill.minimize_length())
}

#[export]
fn zinside(s: f64, delta: f64) -> NumericArray<f64> {
    //let z = |s: &f64| [*s, (PI * *s).sin()];
    //let z = |s: &f64|  [(1.25 + s.cos()) * s.cos(), (1.25 + s.cos()) * s.sin()];
    let z = |s: &f64|  [2.0 * s.cos(), 0.5 * s.sin()];
    //let z = |s: &f64|  [s.cos(), s.sin()];
    let bill = Billiard::new(0.0, 2.0 * PI, 1, delta, z);
    array_to_na(bill.z_inside(&s))
}

#[export]
fn zoutside(s: f64, delta: f64) -> NumericArray<f64> {
    //let z = |s: &f64| [*s, (PI * *s).sin()];
    //let z = |s: &f64|  [-1.25*s.cos() - (2.0 * s).cos(), (-1.25 - 2.0*s.cos()) * s.sin()];
    //let z = |s: &f64|  [(1.25 + s.cos()) * s.cos(), (1.25 + s.cos()) * s.sin()];
    let z = |s: &f64|  [2.0 * s.cos(), 0.5 * s.sin()];
    let bill = Billiard::new(0.0, 2.0 * PI, 1, delta, z);
    array_to_na(bill.z_outside(&s))
}

fn vec_to_na(v: Vec<[f64;2]>) -> NumericArray<f64> {
    NumericArray::from_array(
        &[ v.len(), 2],
        v.into_flattened().as_slice()
    )
}

fn array_to_na(a: [f64;2]) -> NumericArray<f64> {
    NumericArray::from_array(
        &[1, 2],
        &a
    )
}
