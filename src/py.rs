use std::collections::HashMap;

use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};

use crate::{
    collect_stats, order::T2, schemes::greedy::GreedyP, Alternating, AntiLex, Lex, RandomO,
    ThresholdABB, ABB,
};

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
        // Lex orders
        "Random" => Box::new(schemes::RM(())),
        "Lex" => Box::new(schemes::M(Lex)),
        "Alternating" => Box::new(schemes::M(Alternating)),
        "ABB" => Box::new(schemes::M(ABB)),
        "ABB2" => Box::new(schemes::M((ABB, Lex))),
        "ThresholdABB" => Box::new(schemes::M(ThresholdABB {
            thr: get(args, "thr")? as u8,
        })),
        "ThresholdABB2" => Box::new(schemes::M((
            ThresholdABB {
                thr: get(args, "thr")? as u8,
            },
            AntiLex,
        ))),
        "T2" => Box::new(schemes::M((
            T2 {
                thr: get(args, "thr")? as u8,
            },
            AntiLex,
        ))),
        "Thresholds" => Box::new(schemes::M((
            (
                (
                    ThresholdABB {
                        thr: get(args, "thr")? as u8,
                    },
                    ThresholdABB {
                        thr: 2 * get(args, "thr")? as u8,
                    },
                ),
                ThresholdABB {
                    thr: 3 * get(args, "thr")? as u8,
                },
            ),
            AntiLex,
        ))),
        "AntiLex" => Box::new(schemes::M(AntiLex)),
        "Greedy" => Box::new(schemes::M(GreedyP)),

        // Decycling
        "Decycling" => Box::new(schemes::RM(Decycling { double: false })),
        "DoubleDecycling" => Box::new(schemes::RM(Decycling { double: true })),

        // Selection schemes
        "BdAnchor" => Box::new(schemes::BdAnchor { r: get(args, "r")? }),

        // SUS
        "SusLex" => Box::new(schemes::SusAnchor(Lex)),
        "SusAlternating" => Box::new(schemes::SusAnchor(Alternating)),
        "SusABB" => Box::new(schemes::SusAnchor(ABB)),
        "SusABB2" => Box::new(schemes::SusAnchor((ABB, Lex))),
        "SusThresholdABB" => Box::new(schemes::SusAnchor(ThresholdABB {
            thr: get(args, "thr")? as u8,
        })),
        "SusThresholdABB2" => Box::new(schemes::SusAnchor((
            ThresholdABB {
                thr: get(args, "thr")? as u8,
            },
            Lex,
        ))),
        "SusAntiLex" => Box::new(schemes::SusAnchor(AntiLex)),

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
                    miniception_r: false,
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
                    miniception_r: false,
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
) -> PyResult<(f64, Vec<f64>, Vec<f64>, Vec<Vec<f64>>, HashMap<u64, usize>)> {
    let scheme = get_scheme(tp, params)?.build(w, k, sigma);
    let stats = collect_stats(w, k, &text, &*scheme);
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
    stats(tp, text, w, k, sigma, params).map(|(density, _, _, _, _)| density)
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
