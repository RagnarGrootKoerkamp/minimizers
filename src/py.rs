use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};

use crate::{collect_stats, AntiLex, RandomO};

fn get(dict: Option<&Bound<'_, PyDict>>, key: &str) -> PyResult<usize> {
    Ok(dict
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "Missing minimizer parameter {key}. Add {key}=<val>."
            ))
        })?
        .get_item(key)?
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "Missing minimizer parameter {key}. Add {key}=<val>."
            ))
        })?
        .extract()?)
}

fn get_bool(dict: Option<&Bound<'_, PyDict>>, key: &str) -> bool {
    get(dict, key).map_or(false, |x| x == 1)
}

fn get_scheme(tp: &str, args: Option<&Bound<'_, PyDict>>) -> PyResult<Box<dyn super::Params>> {
    use super::schemes;
    use super::schemes::*;
    let mut params: Box<dyn super::Params> = match tp {
        "Rot" => Box::new(RM(RotMinimizer)),
        "AltRot" => Box::new(RM(AltRotMinimizer)),
        "Bruteforce" => Box::new(schemes::BruteforceP),
        "Random" => Box::new(schemes::RM(())),
        "AntiLex" => Box::new(schemes::M(AntiLex)),
        "Decycling" => Box::new(schemes::RM(Decycling { double: false })),
        "DoubleDecycling" => Box::new(schemes::RM(Decycling { double: true })),
        "BdAnchor" => Box::new(schemes::BdAnchor { r: get(args, "r")? }),
        "SusAnchorLex" => Box::new(schemes::SusAnchorLex),
        "SusAnchorALex" => Box::new(schemes::SusAnchorALex),
        "FracMin" => Box::new(schemes::RM(schemes::FracMin { f: get(args, "f")? })),
        "OpenClosed" => {
            if get_bool(args, "anti_lex") {
                Box::new(RM(OpenClosed {
                    r: get(args, "r")?,
                    open: get_bool(args, "open"),
                    closed: get_bool(args, "closed"),
                    open_by_tmer: get_bool(args, "open_tmer"),
                    closed_by_tmer: get_bool(args, "closed_tmer"),
                    other_by_tmer: get_bool(args, "other_tmer"),
                    offset: get(args, "offset").ok(),
                    modulo: get_bool(args, "modulo"),
                    anti_tmer: get_bool(args, "anti_tmer"),
                    o: AntiLex,
                }))
            } else {
                Box::new(RM(OpenClosed {
                    r: get(args, "r")?,
                    open: get_bool(args, "open"),
                    closed: get_bool(args, "closed"),
                    open_by_tmer: get_bool(args, "open_tmer"),
                    closed_by_tmer: get_bool(args, "closed_tmer"),
                    other_by_tmer: get_bool(args, "other_tmer"),
                    offset: get(args, "offset").ok(),
                    modulo: get_bool(args, "modulo"),
                    anti_tmer: get_bool(args, "anti_tmer"),
                    o: RandomO,
                }))
            }
        }
        "Threshold" => Box::new(schemes::ThresholdMinimizerP {
            t: get(args, "t")?,
            h: get(args, "h")?,
            loose: get_bool(args, "loose"),
            open: get_bool(args, "open"),
        }),
        _ => PyResult::Err(PyValueError::new_err("Invalid minimizer type"))?,
    };
    if get_bool(args, "mod") {
        params = Box::new(schemes::ModP {
            r: get(args, "r").unwrap_or(1),
            lr: get_bool(args, "lr"),
            t: get(args, "sampling").unwrap_or(0),
            params,
        });
    }
    Ok(params)
}

#[pyfunction]
#[pyo3(signature = (tp, text, w, k, sigma, **params))]
fn stats(
    tp: &str,
    text: Vec<u8>,
    w: usize,
    k: usize,
    sigma: usize,
    params: Option<&Bound<'_, PyDict>>,
) -> PyResult<(f64, Vec<f64>, Vec<f64>, Vec<Vec<f64>>)> {
    let scheme = get_scheme(tp, params)?.build(w, k, sigma);
    let stats = collect_stats(w, &text, &*scheme);
    Ok(stats)
}

#[pyfunction]
#[pyo3(signature = (tp, text, w, k, l, sigma, **params))]
fn cycle_stats(
    tp: &str,
    text: Vec<u8>,
    w: usize,
    k: usize,
    l: usize,
    sigma: usize,
    params: Option<&Bound<'_, PyDict>>,
) -> PyResult<(f64, Vec<f64>)> {
    let scheme = get_scheme(tp, params)?.build(w, k, sigma);
    let stats = super::cycle_stats(l, &text, &*scheme);
    Ok(stats)
}

#[pyfunction]
#[pyo3(signature = (tp, text, w, k, sigma, **params))]
fn density(
    tp: &str,
    text: Vec<u8>,
    w: usize,
    k: usize,
    sigma: usize,
    params: Option<&Bound<'_, PyDict>>,
) -> PyResult<f64> {
    stats(tp, text, w, k, sigma, params).map(|(density, _, _, _)| density)
}

#[pyfunction]
pub fn generate_random_string(n: usize, sigma: usize) -> PyResult<Vec<u8>> {
    Ok(super::generate_random_string(n, sigma))
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn minimizers(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(stats, m)?)?;
    m.add_function(wrap_pyfunction!(cycle_stats, m)?)?;
    m.add_function(wrap_pyfunction!(density, m)?)?;
    m.add_function(wrap_pyfunction!(generate_random_string, m)?)?;
    Ok(())
}
