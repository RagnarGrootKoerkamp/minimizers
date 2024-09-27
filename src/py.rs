use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};

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
    let scheme: super::MinimizerType = match tp {
        "LrMinimizer" | "RotMinimizer" | "AltRotMinimizer" | "DecyclingMinimizer"
        | "Bruteforce" => {
            serde_json::from_str(&format!("{{\"minimizer_type\": \"{tp}\"}}")).unwrap()
        }
        "Minimizer" => super::MinimizerType::Minimizer {
            ao: get(params, "ao").map_or(false, |x| x == 1),
        },
        "DoubleDecyclingMinimizer" => super::MinimizerType::DoubleDecyclingMinimizer {
            ao: get(params, "ao").map_or(false, |x| x == 1),
        },
        "ModMinimizer" => super::MinimizerType::ModMinimizer {
            r: get(params, "r")?,
            aot: get(params, "aot").map_or(false, |x| x == 1),
        },
        "BdAnchor" => super::MinimizerType::BdAnchor {
            r: get(params, "r")?,
        },
        "SusAnchor" => super::MinimizerType::SusAnchor {
            ao: get(params, "ao").map_or(false, |x| x == 1),
        },
        "Miniception" => super::MinimizerType::Miniception {
            k0: get(params, "k0")?,
            ao: get(params, "ao").map_or(false, |x| x == 1),
            aot: get(params, "aot").map_or(false, |x| x == 1),
        },
        "MiniceptionNew" => super::MinimizerType::MiniceptionNew {
            k0: get(params, "k0")?,
        },
        "ModSampling" => super::MinimizerType::ModSampling {
            k0: get(params, "k0")?,
        },
        "OpenSyncmerMinimizer" => super::MinimizerType::OpenSyncmerMinimizer {
            t: get(params, "t")?,
        },
        "ClosedSyncmerMinimizer" => super::MinimizerType::ClosedSyncmerMinimizer {
            t: get(params, "t")?,
            h: get(params, "h")?,
            loose: get(params, "loose")? == 1,
            open: get(params, "open")? == 1,
        },
        "OpenClosedSyncmerMinimizer" => super::MinimizerType::OpenClosedSyncmerMinimizer {
            t: get(params, "t")?,
        },
        "FracMin" => super::MinimizerType::FracMin {
            f: get(params, "f")?,
        },
        "OcModMinimizer" => super::MinimizerType::OcModMinimizer {
            t: get(params, "t")?,
            offset: get(params, "offset")?,
            use_closed: get(params, "use_closed")? == 1,
            prefer_prefix: get(params, "prefer_prefix")? == 1,
            open_tmer: get(params, "open_tmer")? == 1,
            closed_tmer: get(params, "closed_tmer")? == 1,
            other_tmer: get(params, "other_tmer")? == 1,
            ao: get(params, "ao").map_or(false, |x| x == 1),
            aot: get(params, "aot").map_or(false, |x| x == 1),
        },
        _ => PyResult::Err(PyValueError::new_err("Invalid minimizer type"))?,
    };
    let stats = scheme.collect_stats(&text, w, k, sigma);
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
    let scheme: super::MinimizerType = match tp {
        "LrMinimizer" | "RotMinimizer" | "AltRotMinimizer" | "DecyclingMinimizer"
        | "Bruteforce" => {
            serde_json::from_str(&format!("{{\"minimizer_type\": \"{tp}\"}}")).unwrap()
        }
        "Minimizer" => super::MinimizerType::Minimizer {
            ao: get(params, "ao").map_or(false, |x| x == 1),
        },
        "DoubleDecyclingMinimizer" => super::MinimizerType::DoubleDecyclingMinimizer {
            ao: get(params, "ao").map_or(false, |x| x == 1),
        },
        "ModMinimizer" => super::MinimizerType::ModMinimizer {
            r: get(params, "r")?,
            aot: get(params, "aot").map_or(false, |x| x == 1),
        },
        "BdAnchor" => super::MinimizerType::BdAnchor {
            r: get(params, "r")?,
        },
        "SusAnchor" => super::MinimizerType::SusAnchor {
            ao: get(params, "ao").map_or(false, |x| x == 1),
            modulo: get(params, "modulo").map_or(false, |x| x == 1),
        },
        "Miniception" => super::MinimizerType::Miniception {
            k0: get(params, "k0")?,
            ao: get(params, "ao").map_or(false, |x| x == 1),
            aot: get(params, "aot").map_or(false, |x| x == 1),
        },
        "MiniceptionNew" => super::MinimizerType::MiniceptionNew {
            k0: get(params, "k0")?,
        },
        "ModSampling" => super::MinimizerType::ModSampling {
            k0: get(params, "k0")?,
        },
        "OpenSyncmerMinimizer" => super::MinimizerType::OpenSyncmerMinimizer {
            t: get(params, "t")?,
        },
        "ClosedSyncmerMinimizer" => super::MinimizerType::ClosedSyncmerMinimizer {
            t: get(params, "t")?,
            h: get(params, "h")?,
            loose: get(params, "loose")? == 1,
            open: get(params, "open")? == 1,
        },
        "OpenClosedSyncmerMinimizer" => super::MinimizerType::OpenClosedSyncmerMinimizer {
            t: get(params, "t")?,
        },
        "FracMin" => super::MinimizerType::FracMin {
            f: get(params, "f")?,
        },
        "OcModMinimizer" => super::MinimizerType::OcModMinimizer {
            t: get(params, "t")?,
            offset: get(params, "offset")?,
            use_closed: get(params, "use_closed")? == 1,
            prefer_prefix: get(params, "prefer_prefix")? == 1,
            open_tmer: get(params, "open_tmer")? == 1,
            closed_tmer: get(params, "closed_tmer")? == 1,
            other_tmer: get(params, "other_tmer")? == 1,
            ao: get(params, "ao").map_or(false, |x| x == 1),
            aot: get(params, "aot").map_or(false, |x| x == 1),
        },
        _ => PyResult::Err(PyValueError::new_err("Invalid minimizer type"))?,
    };
    let stats = scheme.cycle_stats(&text, w, k, l, sigma);
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

#[pyfunction]
#[pyo3(signature = (w, k, l))]
fn cycle_partition_lp(w: usize, k: usize, l: usize) -> PyResult<f64> {
    Ok(super::cycle_partition_lp::cycle_partition_lp(
        w,
        k,
        l..l + 1,
    ))
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
    m.add_function(wrap_pyfunction!(cycle_partition_lp, m)?)?;
    Ok(())
}
