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
fn density(
    tp: &str,
    text: Vec<u8>,
    w: usize,
    k: usize,
    sigma: usize,
    params: Option<&Bound<'_, PyDict>>,
) -> PyResult<f64> {
    let scheme: super::MinimizerType = match tp {
        "Minimizer"
        | "LrMinimizer"
        | "ModMinimizer"
        | "RotMinimizer"
        | "AltRotMinimizer"
        | "DecyclingMinimizer"
        | "DoubleDecyclingMinimizer"
        | "Bruteforce" => {
            serde_json::from_str(&format!("{{\"minimizer_type\": \"{tp}\"}}")).unwrap()
        }
        "BdAnchor" => super::MinimizerType::BdAnchor {
            r: get(params, "r")?,
        },
        "Miniception" => super::MinimizerType::Miniception {
            k0: get(params, "k0")?,
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
        "OpenClosedSyncmerMinimizer" => super::MinimizerType::OpenClosedSyncmerMinimizer {
            t: get(params, "t")?,
        },
        "FracMin" => super::MinimizerType::FracMin {
            f: get(params, "f")?,
        },
        _ => PyResult::Err(PyValueError::new_err("Invalid minimizer type"))?,
    };
    let density = scheme.stats(&text, w, k, sigma).0;
    Ok(density)
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
    m.add_function(wrap_pyfunction!(density, m)?)?;
    m.add_function(wrap_pyfunction!(generate_random_string, m)?)?;
    Ok(())
}
