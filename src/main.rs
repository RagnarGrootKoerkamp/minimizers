#![feature(type_alias_impl_trait)]
use std::{
    io::Write,
    path::PathBuf,
    sync::{atomic::AtomicUsize, Mutex},
};

use clap::Parser;
use itertools::Itertools;
use minimizers::{de_bruijn_seq::de_bruijn_sequence, *};
use order::RandomOrder;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde_derive::Serialize;

/// Generate a random string.
fn generate_random_string(n: usize, sigma: usize) -> Vec<u8> {
    (0..n)
        .map(|_| (rand::random::<usize>() % sigma) as u8)
        .collect()
}

/// Returns:
/// - density
/// - position distribution
/// - distance distribution
/// - transfer distribution
fn collect_stats(
    w: usize,
    text: &[u8],
    scheme: impl SamplingScheme,
) -> (f64, Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let it = scheme.stream(text);

    let mut n = 0;
    let mut anchors = 0;
    let mut ps = vec![0; w];
    let mut ds = vec![0; 2 * w + 1];
    let mut transfer = vec![vec![0; w]; w];
    let mut last = 0;
    for (i, idx) in it.enumerate() {
        assert!(
            i <= idx && idx < i + w,
            "Sampled index not in range: {i}<={idx}<{}",
            i + w
        );
        n += 1;
        ps[idx - i] += 1;
        transfer[last - (i - 1)][idx - i] += 1;
        if idx != last {
            anchors += 1;
            ds[w + idx - last] += 1;
            last = idx;
        }
    }
    let density = anchors as f64 / n as f64;
    let ps = ps.into_iter().map(|c| (c * w) as f64 / n as f64).collect();
    let ds = ds
        .into_iter()
        .map(|d| (d * w) as f64 / anchors as f64)
        .collect();
    let transfer = transfer
        .into_iter()
        .map(|row| row.into_iter().map(|c| c as f64 / n as f64).collect())
        .collect();
    (density, ps, ds, transfer)
}

#[derive(Clone, Copy, clap::Subcommand, Debug, Serialize, Default)]
#[serde(tag = "minimizer_type")]
enum MinimizerType {
    #[default]
    Minimizer,
    BdAnchor {
        r: usize,
    },
    Miniception {
        k0: usize,
    },
    MiniceptionNew {
        k0: usize,
    },
    ModSampling {
        k0: usize,
    },
    LrMinimizer,
    ModMinimizer,
    RotMinimizer,
    AltRotMinimizer,
    DecyclingMinimizer,
    DoubleDecyclingMinimizer,
    Bruteforce,
    OpenSyncmerMinimizer {
        t: usize,
        tiebreak: bool,
        transfer: usize,
    },
    OpenClosedSyncmerMinimizer {
        t: usize,
        tiebreak: bool,
        transfer: usize,
    },
    FracMin {
        f: usize,
    },
}

#[derive(Clone, Debug, Serialize, Default)]
struct Result {
    sigma: usize,
    k: usize,
    w: usize,
    l: usize,
    tp: MinimizerType,
    density: f64,
    positions: Vec<f64>,
    dists: Vec<f64>,
    transfer: Vec<Vec<f64>>,
}

/// TODO: Analyze non-forward schemes.
impl MinimizerType {
    #[inline(never)]
    fn stats(
        &self,
        text: &[u8],
        w: usize,
        k: usize,
        sigma: usize,
    ) -> (f64, Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
        let o = RandomOrder;
        match self {
            MinimizerType::Minimizer => collect_stats(w, text, Minimizer::new(k, w, o)),
            MinimizerType::BdAnchor { r } => collect_stats(w, text, BdAnchor::new(w + *r, *r)),
            MinimizerType::Miniception { k0 } => {
                collect_stats(w, text, Miniception::new(w, k, *k0, o))
            }
            MinimizerType::MiniceptionNew { k0 } => {
                collect_stats(w, text, MiniceptionNew::new(w, k, *k0, o))
            }
            MinimizerType::ModSampling { k0 } => {
                collect_stats(w, text, ModSampling::new(k, w, *k0, o))
            }
            MinimizerType::LrMinimizer => collect_stats(w, text, ModSampling::lr_minimizer(k, w)),
            MinimizerType::ModMinimizer => collect_stats(w, text, ModSampling::mod_minimizer(k, w)),
            MinimizerType::RotMinimizer => collect_stats(w, text, RotMinimizer::new(k, w, o)),
            MinimizerType::AltRotMinimizer => collect_stats(w, text, AltRotMinimizer::new(k, w, o)),
            MinimizerType::DecyclingMinimizer => {
                collect_stats(w, text, Decycling::new(k, w, o, false))
            }
            MinimizerType::DoubleDecyclingMinimizer => {
                collect_stats(w, text, Decycling::new(k, w, o, true))
            }
            MinimizerType::Bruteforce => {
                let m = bruteforce::bruteforce_minimizer(k, w, sigma).1;
                collect_stats(w, text, m)
            }
            MinimizerType::OpenSyncmerMinimizer {
                t,
                tiebreak,
                transfer,
            } => collect_stats(
                w,
                text,
                OpenSyncmer::new(k, w, *t, *tiebreak, false, *transfer),
            ),
            MinimizerType::OpenClosedSyncmerMinimizer {
                t,
                tiebreak,
                transfer,
            } => collect_stats(
                w,
                text,
                OpenSyncmer::new(k, w, *t, *tiebreak, true, *transfer),
            ),
            MinimizerType::FracMin { f } => collect_stats(w, text, FracMin::new(k, w, *f)),
        }
    }

    fn try_params(&self, w: usize, k: usize) -> Vec<Self> {
        let l = w + k - 1;
        match self {
            MinimizerType::Minimizer
            | MinimizerType::DecyclingMinimizer
            | MinimizerType::DoubleDecyclingMinimizer
            | MinimizerType::Bruteforce => {
                vec![*self]
            }
            MinimizerType::BdAnchor { .. } => {
                let r_max = k;
                (0.min(r_max)..=10.min(r_max))
                    .map(|r| MinimizerType::BdAnchor { r })
                    .collect()
            }
            MinimizerType::Miniception { .. } => {
                let k0_min = (2 * k).saturating_sub(l + 1).max(1);
                let k0_max = k;
                if k0_min > k0_max {
                    return vec![];
                }

                return vec![MinimizerType::Miniception { k0: k0_min }];
                // FIXME
                // let start = 1.max(k0_min);
                // let end = 10.max(k0_min + 2).min(k0_max);
                // (start..=end)
                //     .map(|k0| MinimizerType::Miniception { k0 })
                //     .collect()
            }
            MinimizerType::MiniceptionNew { .. } => {
                let k0_min = (2 * k).saturating_sub(l + 1).max(1);
                let k0_max = k;
                if k0_min > k0_max {
                    return vec![];
                }
                return vec![MinimizerType::MiniceptionNew { k0: k0_min }];

                // FIXME
                // let start = 1.max(k0_min);
                // let end = 10.max(k0_min + 2).min(k0_max);
                // (start..=end)
                //     .map(|k0| MinimizerType::MiniceptionNew { k0 })
                //     .collect()
            }
            MinimizerType::ModSampling { .. } => {
                let k0_min = 1;
                let k0_max = l;
                (k0_min..=k0_max)
                    .map(|k0| MinimizerType::ModSampling { k0 })
                    .collect()
            }
            MinimizerType::LrMinimizer => {
                if k > w + 4 {
                    vec![*self]
                } else {
                    vec![]
                }
            }
            MinimizerType::ModMinimizer => {
                if k > 1 {
                    vec![*self]
                } else {
                    vec![]
                }
            }
            MinimizerType::RotMinimizer => {
                if k % w == 0 {
                    vec![*self]
                } else {
                    vec![]
                }
            }
            MinimizerType::AltRotMinimizer => {
                if k > w {
                    vec![*self]
                } else {
                    vec![]
                }
            }
            MinimizerType::OpenSyncmerMinimizer { .. } => {
                let t_min = 1;
                // FIXME: For large k, t>1 is better. But we focus on small k for now.
                let t_max = k.min(1);
                (t_min..=t_max)
                    .step_by(3)
                    .cartesian_product([true])
                    // .cartesian_product(0..w)
                    .map(|(t, tiebreak)| MinimizerType::OpenSyncmerMinimizer {
                        t,
                        tiebreak,
                        transfer: 0,
                    })
                    .collect()
            }
            MinimizerType::OpenClosedSyncmerMinimizer { .. } => {
                vec![MinimizerType::OpenClosedSyncmerMinimizer {
                    t: 1,
                    tiebreak: true,
                    transfer: 0,
                }]
            }
            MinimizerType::FracMin { .. } => (1..w).map(|f| MinimizerType::FracMin { f }).collect(),
        }
    }
}

#[derive(clap::Subcommand)]
enum Command {
    Run {
        #[clap(short, default_value_t = 40)]
        w: usize,
        #[arg(short, default_value_t = 10)]
        k: usize,
        #[clap(subcommand)]
        tp: MinimizerType,
    },
    Compare {
        #[clap(short, default_value_t = 40)]
        w: usize,
        #[arg(short, default_value_t = 10)]
        k: usize,
    },
    Eval {
        /// Write to file.
        #[clap(short, long)]
        output: Option<PathBuf>,
        /// Run in statistics mode.
        #[clap(short, long)]
        stats: bool,
        #[clap(short, long)]
        small: bool,
    },
}

/// Print expected and measured densities for different w and r.
#[derive(clap::Parser)]
struct Args {
    /// Length of the generated random string.
    #[clap(short, default_value_t = 100000)]
    n: usize,
    /// Alphabet size.
    #[clap(short, long, default_value_t = 4)]
    sigma: usize,
    #[clap(subcommand)]
    command: Command,
}

fn main() {
    let args = Args::parse();
    let text = &generate_random_string(args.n, args.sigma);

    match args.command {
        Command::Run { tp, w, k } => {
            eprintln!("Running {tp:?}:");
            let (d, _ps, _ds, _transfer) = tp.stats(text, w, k, args.sigma);
            eprintln!("  Density : {d:.3}");
            // eprintln!("  Poss    : {ps:.5?}");
            // eprintln!("  Dists<0 : {:.5?}", &ds[0..w]);
            // eprintln!("  Dists>0 : {:.5?}", &ds[w + 1..]);
            // eprintln!("  Transfer:");
            // for row in &transfer {
            // eprintln!("    {row:.5?}");
            // }
        }
        Command::Compare { w, k } => compare_methods(text, k, w, args.sigma),

        Command::Eval {
            output,
            stats,
            small,
        } => {
            let mut base_types = vec![
                MinimizerType::Minimizer,
                // MinimizerType::BdAnchor { r: 0 },
                MinimizerType::Miniception { k0: 0 },
                // MinimizerType::MiniceptionNew { k0: 0 },
                // MinimizerType::ModSampling { k0: 0 },
                MinimizerType::LrMinimizer,
                MinimizerType::ModMinimizer,
                MinimizerType::RotMinimizer,
                MinimizerType::AltRotMinimizer,
                MinimizerType::DecyclingMinimizer,
                MinimizerType::DoubleDecyclingMinimizer,
                MinimizerType::OpenSyncmerMinimizer {
                    t: 0,
                    tiebreak: false,
                    transfer: 0,
                },
                // MinimizerType::OpenClosedSyncmerMinimizer {
                //     t: 0,
                //     tiebreak: false,
                //     transfer: 0,
                // },
                // MinimizerType::FracMin { f: 0 },
            ];
            if small {
                base_types.push(MinimizerType::Bruteforce);
            }

            let ks = if small {
                &[1, 2, 3][..]
            } else if stats {
                &[4, 8, 16, 32][..]
            } else {
                &[
                    5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26,
                    28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84,
                    88, 92, 96, 100,
                ][..]
            };
            let ws = if small {
                &[2][..]
            } else if stats {
                &[8, 16, 32, 64][..]
            } else {
                &[8, 24][..]
            };
            let k_w_tp = ks
                .iter()
                .cartesian_product(ws)
                .cartesian_product(base_types)
                .collect_vec();
            let results = Mutex::new(vec![None; k_w_tp.len()]);
            let done = AtomicUsize::new(0);
            let total = k_w_tp.len();
            let text = if small {
                &de_bruijn_sequence(args.sigma, 8)
            } else {
                text
            };
            k_w_tp.par_iter().for_each(|&((&k, &w), tp)| {
                if small && k == 3 && args.sigma == 3 {
                    return;
                }
                let l = w + k - 1;
                let tps = tp.try_params(w, k);
                for tp in tps {
                    let (density, positions, dists, transfer) = tp.stats(text, w, k, args.sigma);
                    let mut results = results.lock().unwrap();
                    results.push(Some(Result {
                        sigma: args.sigma,
                        k,
                        w,
                        l,
                        tp,
                        density,
                        positions,
                        dists,
                        transfer,
                    }));
                    // let done = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    // eprint!("{done:>3}/{total:>3}: k={k} w={w} l={l} tp={tp:?} d={density:.3}\r");
                }
                // let Some(((density, positions, dists, transfer), tp)) = tps
                //     .iter()
                //     .map(|tp| (tp.stats(text, w, k, args.sigma), tp))
                //     .min_by(|(ld, _), (rd, _)| {
                //         if ld.0 < rd.0 {
                //             std::cmp::Ordering::Less
                //         } else {
                //             std::cmp::Ordering::Greater
                //         }
                //     })
                // else {
                //     done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                //     return;
                // };
                // let mut results = results.lock().unwrap();
                // results.push(Some(Result {
                //     sigma: args.sigma,
                //     k,
                //     w,
                //     l,
                //     tp: *tp,
                //     density,
                //     positions,
                //     dists,
                //     transfer,
                // }));
                let done = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                eprint!("{done:>3}/{total:>3}: k={k} w={w} l={l} tp={tp:?}\r");
            });
            eprintln!();

            if let Some(output) = output {
                let result_json = serde_json::to_string(
                    &results
                        .into_inner()
                        .unwrap()
                        .into_iter()
                        .filter_map(|x| x)
                        .collect_vec(),
                )
                .unwrap();
                let mut file = std::fs::File::create(output).unwrap();
                file.write_all(result_json.as_bytes()).unwrap();
            }
        }
    }
}

fn compare_methods(text: &[u8], k: usize, w: usize, sigma: usize) {
    let base_types = vec![
        MinimizerType::Minimizer,
        // MinimizerType::BdAnchor { r: 0 },
        MinimizerType::Miniception { k0: 0 },
        MinimizerType::MiniceptionNew { k0: 0 },
        MinimizerType::ModSampling { k0: 0 },
        MinimizerType::RotMinimizer,
        MinimizerType::DecyclingMinimizer,
        MinimizerType::DoubleDecyclingMinimizer,
    ];

    base_types.par_iter().for_each(|&tp| {
        let l = w + k - 1;
        let tps = tp.try_params(w, k);
        let Some(((density, _positions, _dists, _transfer), tp)) = tps
            .iter()
            .map(|tp| (tp.stats(text, w, k, sigma), tp))
            .min_by(|(ld, _), (rd, _)| {
                if ld.0 < rd.0 {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Greater
                }
            })
        else {
            return;
        };
        eprintln!("k={k} w={w} l={l} d={density:.3} tp={tp:?}");
    });

    let best = bruteforce::bruteforce_minimizer(k, w, sigma);
    let d = (best.0 .0) as f32 / best.0 .1 as f32;
    eprintln!("k={k} w={w} l={} d={d:.3} tp=Bruteforce", k + w - 1);
}

#[cfg(test)]
mod test {

    use super::*;
    #[test]
    fn minimizers() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20 {
            for w in 1..=20 {
                let m = RandomMinimizer::new(k, w, RandomOrder);
                let stream = m.stream(&text).collect_vec();
                let stream_naive = m.stream_naive(&text).collect_vec();
                assert_eq!(stream, stream_naive);
            }
        }
    }

    #[test]
    fn mod_minimizers() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20 {
            for w in 1..=20 {
                let l = k + w - 1;
                for t in 1..=l {
                    let m = ModSampling::new(k, w, t, RandomOrder);
                    let stream = m.stream(&text).collect_vec();
                    let stream_naive = m.stream_naive(&text).collect_vec();
                    assert_eq!(stream, stream_naive);
                }
            }
        }
    }

    #[test]
    fn miniception() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20usize {
            for w in 1..=20 {
                for k0 in k.saturating_sub(w).max(1)..=k {
                    let m = Miniception::new(w, k, k0, RandomOrder);
                    let stream = m.stream(&text).collect_vec();
                    let stream_naive = m.stream_naive(&text).collect_vec();
                    assert_eq!(stream, stream_naive);
                }
            }
        }
    }

    #[test]
    fn miniception_new() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20usize {
            for w in 1..=20 {
                for k0 in k.saturating_sub(w).max(1)..=k {
                    let m = MiniceptionNew::new(w, k, k0, RandomOrder);
                    let stream = m.stream(&text).collect_vec();
                    let stream_naive = m.stream_naive(&text).collect_vec();
                    assert_eq!(stream, stream_naive);
                }
            }
        }
    }

    #[test]
    fn syncmer_minimizer() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20 {
            for w in 1..=20 {
                for t in 1..k {
                    eprintln!("k {k} w {w} t {t}");
                    let m = OpenSyncmer::new(k, w, t, true, false, 0);
                    let stream = m.stream(&text).collect_vec();
                    let stream_naive = m.stream_naive(&text).collect_vec();
                    assert_eq!(stream, stream_naive);
                }
            }
        }
    }
}
