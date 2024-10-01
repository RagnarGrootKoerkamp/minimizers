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
        "RotMinimizer" => Box::new(schemes::RotMinimizerP),
        "AltRotMinimizer" => Box::new(schemes::AltRotMinimizerP),
        "Bruteforce" => Box::new(schemes::BruteforceP),
        "Random" => Box::new(schemes::M(RandomO)),
        "AntiLex" => Box::new(schemes::M(AntiLex)),
        "DecyclingMinimizer" => Box::new(schemes::M((Decycling { double: false }, RandomO))),
        "DoubleDecyclingMinimizer" => Box::new(schemes::M((Decycling { double: true }, RandomO))),
        "BdAnchor" => Box::new(schemes::BdAnchorP { r: get(args, "r")? }),
        "SusAnchor" => Box::new(schemes::SusAnchorP {
            ao: get_bool(args, "ao"),
        }),
        // TODO: variants for the two order arguments.
        "Miniception" => Box::new(M((
            Miniception {
                r: get(args, "k0")?,
                o: RandomO,
            },
            RandomO,
        ))),
        "MiniceptionNew" => Box::new(schemes::MiniceptionNewP {
            k0: get(args, "k0")?,
        }),
        "OpenSyncmerMinimizer" => Box::new(schemes::OpenSyncmerMinimizerP { t: get(args, "t")? }),
        "ClosedSyncmerMinimizer" => Box::new(schemes::ThresholdMinimizerP {
            t: get(args, "t")?,
            h: get(args, "h")?,
            loose: get_bool(args, "loose"),
            open: get_bool(args, "open"),
        }),
        "OpenClosedSyncmerMinimizer" => {
            Box::new(schemes::OpenClosedSyncmerMinimizerP { t: get(args, "t")? })
        }
        "FracMin" => Box::new(schemes::M(schemes::FracMin { f: get(args, "f")? })),
        "OcModMinimizer" => Box::new(schemes::OcModMinimizerP {
            t: get(args, "t")?,
            offset: get(args, "offset")?,
            use_closed: get_bool(args, "use_closed"),
            prefer_prefix: get_bool(args, "prefer_prefix"),
            open_tmer: get_bool(args, "open_tmer"),
            closed_tmer: get_bool(args, "closed_tmer"),
            other_tmer: get_bool(args, "other_tmer"),
            ao: get_bool(args, "ao"),
            aot: get_bool(args, "aot"),
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
