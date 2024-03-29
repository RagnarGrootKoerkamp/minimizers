#![feature(exclusive_range_pattern, type_alias_impl_trait)]
use std::{io::Write, path::PathBuf, sync::atomic::AtomicUsize};

use clap::Parser;
use itertools::Itertools;
use minimizers::*;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde_derive::Serialize;

/// Generate a random string.
fn generate_random_string(n: usize, sigma: usize) -> Vec<u8> {
    (0..n)
        .map(|_| (rand::random::<usize>() % sigma) as u8)
        .collect()
}

/// scheme must return a value in [0, l - k].
///
/// Returns anchors and position distribution.
fn stream<'a>(
    text: &'a [u8],
    w: usize,
    k: usize,
    mut scheme: impl FnMut(&[u8]) -> usize + 'a,
) -> impl Iterator<Item = usize> + 'a {
    text.windows(w + k - 1)
        .enumerate()
        .map(move |(i, w)| i + scheme(w))
}

/// Returns:
/// - density
/// - position distribution
/// - distance distribution
/// - transfer distribution
fn collect_stats(
    w: usize,
    it: impl Iterator<Item = usize>,
) -> (f64, Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let mut n = 0;
    let mut anchors = 0;
    let mut ps = vec![0; w];
    let mut ds = vec![0; 2 * w + 1];
    let mut transfer = vec![vec![0; w]; w];
    let mut last = 0;
    for (i, idx) in it.enumerate() {
        assert!(i <= idx && idx < i + w);
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
    BiMinimizer,
    ModMinimizer {
        k0: usize,
    },
    LrMinimizer {
        k0: usize,
    },
    RotMinimizer,
    DecyclingMinimizer,
    DoubleDecyclingMinimizer,
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
    fn stats(&self, text: &[u8], w: usize, k: usize) -> (f64, Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
        match self {
            MinimizerType::Minimizer => collect_stats(w, text_random_minimizers(text, w, k)),
            MinimizerType::BdAnchor { r } => {
                collect_stats(w, stream(text, w, k, |lmer| bd_anchor(lmer, *r)))
            }
            MinimizerType::Miniception { k0 } => {
                collect_stats(w, text_miniception(text, w, k, *k0))
            }
            MinimizerType::MiniceptionNew { k0 } => {
                collect_stats(w, text_miniception_new(text, w, k, *k0))
            }
            MinimizerType::BiMinimizer => {
                let last = &mut 0;
                collect_stats(
                    w,
                    stream(text, w, k, move |lmer| robust_biminimizer(lmer, k, last)),
                )
            }
            MinimizerType::ModMinimizer { k0 } => {
                collect_stats(w, text_mod_minimizers(text, w, k, *k0))
            }
            MinimizerType::LrMinimizer { k0 } => {
                collect_stats(w, text_lr_minimizers(text, w, k, *k0))
            }
            MinimizerType::RotMinimizer => {
                collect_stats(w, stream(text, w, k, move |lmer| rot_minimizer(lmer, k)))
            }
            MinimizerType::DecyclingMinimizer => {
                let cs = decycling_minimizer_init(k);
                collect_stats(
                    w,
                    stream(text, w, k, move |lmer| decycling_minimizer(lmer, k, &cs)),
                )
            }
            MinimizerType::DoubleDecyclingMinimizer => {
                let cs = decycling_minimizer_init(k);
                collect_stats(
                    w,
                    stream(text, w, k, move |lmer| {
                        double_decycling_minimizer(lmer, k, &cs)
                    }),
                )
            }
        }
    }

    fn try_params(&self, w: usize, k: usize) -> Vec<Self> {
        let l = w + k - 1;
        match self {
            MinimizerType::Minimizer
            | MinimizerType::BiMinimizer
            | MinimizerType::DecyclingMinimizer
            | MinimizerType::DoubleDecyclingMinimizer => {
                vec![*self]
            }
            MinimizerType::BdAnchor { .. } => {
                let r_max = k;
                (0.min(r_max)..=10.min(r_max))
                    .map(|r| MinimizerType::BdAnchor { r })
                    .collect()
            }
            MinimizerType::Miniception { .. } => {
                let k0_min = (2 * k).saturating_sub(l + 1);
                let k0_max = k;
                if k0_min > k0_max {
                    return vec![];
                }

                let start = 1.max(k0_min);
                let end = 10.max(k0_min + 2).min(k0_max);
                (start..=end)
                    .map(|k0| MinimizerType::Miniception { k0 })
                    .collect()
            }
            MinimizerType::MiniceptionNew { .. } => {
                let k0_min = (2 * k).saturating_sub(l + 1);
                let k0_max = k;
                if k0_min > k0_max {
                    return vec![];
                }

                let start = 1.max(k0_min);
                let end = 10.max(k0_min + 2).min(k0_max);
                (start..=end)
                    .map(|k0| MinimizerType::MiniceptionNew { k0 })
                    .collect()
            }
            MinimizerType::ModMinimizer { .. } => {
                let k0_min = 1;
                let k0_max = l;
                (k0_min..=k0_max)
                    .map(|k0| MinimizerType::ModMinimizer { k0 })
                    .collect()
            }
            MinimizerType::LrMinimizer { .. } => {
                // k <= (l+k0+1)/2
                // 2k <= l + k0 + 1
                // 2k - l - 1 <= k0
                let k0_min = 1.max((2 * k - 1).saturating_sub(l));
                let k0_max = k;
                (k0_min..=k0_max)
                    .map(|k0| MinimizerType::LrMinimizer { k0 })
                    .collect()
            }
            MinimizerType::RotMinimizer => {
                if k > w {
                    vec![*self]
                } else {
                    vec![]
                }
            }
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
    Eval {
        /// Write to file.
        #[clap(short, long)]
        output: Option<PathBuf>,
        /// Run in statistics mode.
        #[clap(short, long)]
        stats: bool,
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
            let (d, ps, ds, transfer) = tp.stats(text, w, k);
            eprintln!("  Density : {d:.3}");
            eprintln!("  Poss    : {ps:.5?}");
            eprintln!("  Dists<0 : {:.5?}", &ds[0..w]);
            eprintln!("  Dists>0 : {:.5?}", &ds[w + 1..]);
            eprintln!("  Transfer:");
            for row in &transfer {
                eprintln!("    {row:.5?}");
            }
        }
        Command::Eval { output, stats } => {
            let base_types = [
                MinimizerType::Minimizer,
                // MinimizerType::BdAnchor { r: 0 },
                MinimizerType::Miniception { k0: 0 },
                MinimizerType::MiniceptionNew { k0: 0 },
                // MinimizerType::BiMinimizer,
                MinimizerType::LrMinimizer { k0: 0 },
                MinimizerType::ModMinimizer { k0: 0 },
                MinimizerType::RotMinimizer,
                MinimizerType::DecyclingMinimizer,
                MinimizerType::DoubleDecyclingMinimizer,
            ];

            let ks = if stats {
                &[4, 8, 16, 32][..]
            } else {
                &[
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                    23, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 52, 56, 60, 64,
                ][..]
            };
            let ws = if stats {
                &[8, 16, 32, 64][..]
            } else {
                &[8, 16, 32][..]
            };
            let k_w_tp = ks
                .iter()
                .cartesian_product(ws)
                .cartesian_product(base_types)
                .collect_vec();
            let mut results = vec![Result::default(); k_w_tp.len()];
            let done = AtomicUsize::new(0);
            let total = k_w_tp.len();
            k_w_tp
                .par_iter()
                .zip(&mut results)
                .for_each(|(&((&k, &w), tp), result)| {
                    let l = w + k - 1;
                    let tps = tp.try_params(w, k);
                    let Some(((density, positions, dists, transfer), tp)) = tps
                        .iter()
                        .map(|tp| (tp.stats(text, w, k), tp))
                        .min_by(|(ld, _), (rd, _)| {
                            if ld.0 < rd.0 {
                                std::cmp::Ordering::Less
                            } else {
                                std::cmp::Ordering::Greater
                            }
                        })
                    else {
                        done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        return;
                    };
                    *result = Result {
                        sigma: args.sigma,
                        k,
                        w,
                        l,
                        tp: *tp,
                        density,
                        positions,
                        dists,
                        transfer,
                    };
                    let done = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    eprint!("{done:>3}/{total:>3}: k={k} w={w} l={l} tp={tp:?} d={density:.3}\r");
                });
            eprintln!();

            if let Some(output) = output {
                let result_json = serde_json::to_string(&results).unwrap();
                let mut file = std::fs::File::create(output).unwrap();
                file.write_all(result_json.as_bytes()).unwrap();
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn minimizers() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20 {
            for w in 1..=20 {
                let anchors = stream(&text, w, k, |lmer| random_minimizer(lmer, k))
                    .dedup()
                    .collect_vec();
                let minimizers = text_random_minimizers(&text, w, k).dedup().collect_vec();
                assert_eq!(anchors, minimizers);
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
                    let anchors = stream(&text, w, k, |lmer| mod_minimizer(lmer, k, t))
                        .dedup()
                        .collect_vec();
                    let minimizers = text_mod_minimizers(&text, w, k, t).dedup().collect_vec();
                    assert_eq!(anchors, minimizers);
                }
            }
        }
    }

    #[test]
    fn lr_minimizers() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20usize {
            for w in 1..=20 {
                let l = k + w - 1;
                let k0_min = 1.max((2 * k - 1).saturating_sub(l));
                let k0_max = k;
                for t in k0_min..=k0_max {
                    let anchors = stream(&text, w, k, |lmer| lr_minimizer(lmer, k, t))
                        .dedup()
                        .collect_vec();
                    let minimizers = text_lr_minimizers(&text, w, k, t).dedup().collect_vec();
                    assert_eq!(anchors, minimizers);
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
                    let anchors = stream(&text, w, k, |lmer| super::miniception(lmer, k, k0))
                        .dedup()
                        .collect_vec();
                    let minimizers = text_miniception(&text, w, k, k0).dedup().collect_vec();
                    assert_eq!(anchors, minimizers);
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
                    let anchors = stream(&text, w, k, |lmer| super::miniception_new(lmer, k, k0))
                        .dedup()
                        .collect_vec();
                    let minimizers = text_miniception_new(&text, w, k, k0).dedup().collect_vec();
                    assert_eq!(anchors, minimizers);
                }
            }
        }
    }
}
