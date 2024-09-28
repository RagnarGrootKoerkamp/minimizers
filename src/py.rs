use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};

use crate::collect_stats;

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

fn get_scheme(tp: &str, params: Option<&Bound<'_, PyDict>>) -> PyResult<Box<dyn super::Params>> {
    use super::schemes;
    Ok(match tp {
        "LrMinimizer" | "RotMinimizer" | "AltRotMinimizer" | "DecyclingMinimizer"
        | "Bruteforce" => {
            serde_json::from_str(&format!("{{\"minimizer_type\": \"{tp}\"}}")).unwrap()
        }
        "Minimizer" => Box::new(schemes::MinimizerP {
            ao: get(params, "ao").map_or(false, |x| x == 1),
        }),
        "DoubleDecyclingMinimizer" => Box::new(schemes::DoubleDecyclingP {
            ao: get(params, "ao").map_or(false, |x| x == 1),
        }),
        "ModMinimizer" => Box::new(schemes::ModMinimizerP {
            r: get(params, "r")?,
            aot: get(params, "aot").map_or(false, |x| x == 1),
        }),
        "BdAnchor" => Box::new(schemes::BdAnchorP {
            r: get(params, "r")?,
        }),
        "SusAnchor" => Box::new(schemes::SusAnchorP {
            ao: get(params, "ao").map_or(false, |x| x == 1),
            modulo: get(params, "modulo").map_or(false, |x| x == 1),
        }),
        "Miniception" => Box::new(schemes::MiniceptionP {
            k0: get(params, "k0")?,
            ao: get(params, "ao").map_or(false, |x| x == 1),
            aot: get(params, "aot").map_or(false, |x| x == 1),
        }),
        "MiniceptionNew" => Box::new(schemes::MiniceptionNewP {
            k0: get(params, "k0")?,
        }),
        "ModSampling" => Box::new(schemes::ModSamplingP {
            k0: get(params, "k0")?,
        }),
        "OpenSyncmerMinimizer" => Box::new(schemes::OpenSyncmerMinimizerP {
            t: get(params, "t")?,
        }),
        "ClosedSyncmerMinimizer" => Box::new(schemes::ThresholdMinimizerP {
            t: get(params, "t")?,
            h: get(params, "h")?,
            loose: get(params, "loose")? == 1,
            open: get(params, "open")? == 1,
        }),
        "OpenClosedSyncmerMinimizer" => Box::new(schemes::OpenClosedSyncmerMinimizerP {
            t: get(params, "t")?,
        }),
        "FracMin" => Box::new(schemes::FracMinP {
            f: get(params, "f")?,
        }),
        "OcModMinimizer" => Box::new(schemes::OcModMinimizerP {
            t: get(params, "t")?,
            offset: get(params, "offset")?,
            use_closed: get(params, "use_closed")? == 1,
            prefer_prefix: get(params, "prefer_prefix")? == 1,
            open_tmer: get(params, "open_tmer")? == 1,
            closed_tmer: get(params, "closed_tmer")? == 1,
            other_tmer: get(params, "other_tmer")? == 1,
            ao: get(params, "ao").map_or(false, |x| x == 1),
            aot: get(params, "aot").map_or(false, |x| x == 1),
        }),
        _ => PyResult::Err(PyValueError::new_err("Invalid minimizer type"))?,
    })
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
