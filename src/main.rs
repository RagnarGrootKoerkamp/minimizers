#![feature(exclusive_range_pattern, type_alias_impl_trait)]
use std::cmp::Reverse;

use clap::Parser;
use itertools::Itertools;
use serde_derive::Serialize;

/// Generate a random string.
fn generate_random_string(n: usize, sigma: usize) -> Vec<u8> {
    (0..n)
        .map(|_| (rand::random::<usize>() % sigma) as u8)
        .collect()
}

fn h(kmer: &[u8]) -> u64 {
    fxhash::hash64(kmer)
}

/// Rightmost element with minimal hash.
fn minimizer(s: &[u8], k: usize) -> usize {
    assert!(k > 0);
    assert!(k <= s.len());
    s.windows(k)
        .enumerate()
        .min_by_key(|&(i, kmer)| (h(kmer), Reverse(i)))
        .unwrap()
        .0
}

/// Rightmost element with minimal hash, or last if possible.
fn robust_minimizer(s: &[u8], k: usize, last: &mut usize) -> usize {
    assert!(k <= s.len());
    let (i, kmer) = s
        .windows(k)
        .enumerate()
        .min_by_key(|&(i, kmer)| (h(kmer), Reverse(i)))
        .unwrap();
    if *last > 0 && h(&s[*last - 1..*last - 1 + k]) == h(kmer) {
        *last -= 1;
    } else {
        *last = i;
    }
    *last
}

fn bd_anchor(s: &[u8], r: usize) -> usize {
    assert!(r <= s.len());
    let mut best = 0;
    for i in 1..s.len().saturating_sub(r) {
        for j in 0..s.len() {
            if s[(i + j) % s.len()] != s[(best + j) % s.len()] {
                if s[(i + j) % s.len()] < s[(best + j) % s.len()] {
                    best = i;
                }
                break;
            }
        }
    }
    best
}

/// 1. Select kmers whose k0-minimizer is at their start or end.
/// 2. From these, select a minimal one.
// TODO: Same hash or distinct hashes?
fn miniception(s: &[u8], k: usize, k0: usize) -> usize {
    let w = s.len() - k + 1;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    let w0 = k - k0;
    let _l0 = w0 + k0 - 1;

    s.windows(k)
        .enumerate()
        .filter(|(_, kmer)| {
            let i = minimizer(kmer, k0);
            assert!(i <= k - k0);
            i == 0 || i == k - k0
        })
        .min_by_key(|&(i, w)| (h(w), Reverse(i)))
        .unwrap()
        .0
}

fn robust_biminimizer_bot(s: &[u8], k: usize, last: &mut usize) -> usize {
    let mut vals = s
        .windows(k)
        .enumerate()
        .map(|(i, w)| (h(w), Reverse(i)))
        .collect_vec();
    vals.sort();
    let i1 = vals[0].1 .0;
    if vals.len() == 1 {
        *last = i1;
        return i1;
    }
    let i2 = vals[1].1 .0;
    if *last == i1 + 1 || *last == i2 + 1 {
        *last -= 1;
    } else {
        *last = i1.max(i2);
    }
    *last
}

fn robust_biminimizer(s: &[u8], k: usize, last: &mut usize) -> usize {
    let h2 = |x: &[u8]| fxhash::hash64(&(h(x) ^ 3283342332002053234u64));
    let i1 = s
        .windows(k)
        .enumerate()
        .min_by_key(|&(i, w)| (h(w), Reverse(i)))
        .unwrap()
        .0;
    let i2 = s
        .windows(k)
        .enumerate()
        .min_by_key(|&(i, w)| (h2(w), Reverse(i)))
        .unwrap()
        .0;
    if *last == i1 + 1 || *last == i2 + 1 {
        *last -= 1;
    } else {
        *last = i1.max(i2);
    }
    *last
}

fn reduced_minimizer(s: &[u8], k: usize, k0: usize) -> usize {
    let l = s.len();
    let _w = l - k + 1;
    assert!(k0 <= k);
    assert!(k <= (l + k0 + 1) / 2);
    let idx = s
        .windows(k0)
        .enumerate()
        .min_by_key(|&(i, w)| (fxhash::hash64(w), Reverse(i)))
        .unwrap()
        .0;
    let i = if idx + k <= l { idx } else { idx - (k - k0) };
    assert!(i <= l, "l={l} idx={idx} k={k} k0={k0} i={i}");
    assert!(i + k <= l, "l={l} idx={idx} k={k} k0={k0} i={i}");
    i
}

/// Compute the density of the sampling scheme.
/// A function to select the minimal kmer from an lmer.
/// The length of the slice is l.
/// k must be bound by the function.
///
/// Must return a value in [0, l - k].
fn density(text: &[u8], l: usize, mut scheme: impl FnMut(&[u8]) -> usize) -> f64 {
    // Collect positions.
    let mut anchors = text
        .windows(l)
        .enumerate()
        .map(|(i, w)| i + scheme(w))
        .dedup()
        .collect_vec();

    // Dedup anchors to ensure valid results for non-forward schemes (bd-anchors).
    anchors.sort();
    anchors.dedup();

    anchors.len() as f64 / text.windows(l).len() as f64
}

#[derive(Clone, Copy, clap::Subcommand, Debug, Serialize)]
#[serde(tag = "minimizer_type")]
enum MinimizerType {
    Minimizer,
    RobustMinimizer,
    BdAnchor { r: usize },
    Miniception { k0: usize },
    BiMinimizer,
    BiMinimizerBot,
    ReducedMinimizer { k0: usize },
}

#[derive(Clone, Copy, Debug, Serialize)]
struct Result {
    sigma: usize,
    k: usize,
    w: usize,
    l: usize,
    tp: MinimizerType,
    density: f64,
}

impl MinimizerType {
    fn density(&self, text: &[u8], l: usize, k: usize) -> f64 {
        match self {
            MinimizerType::Minimizer => density(text, l, |lmer| minimizer(lmer, k)),
            MinimizerType::RobustMinimizer => {
                let last = &mut 0;
                density(text, l, |lmer| robust_minimizer(lmer, k, last))
            }
            MinimizerType::BdAnchor { r } => density(text, l, |lmer| bd_anchor(lmer, *r)),
            MinimizerType::Miniception { k0 } => density(text, l, |lmer| miniception(lmer, k, *k0)),
            MinimizerType::BiMinimizer => {
                let last = &mut 0;
                density(text, l, move |lmer| robust_biminimizer(lmer, k, last))
            }
            MinimizerType::BiMinimizerBot => {
                let last = &mut 0;
                density(text, l, move |lmer| robust_biminimizer_bot(lmer, k, last))
            }
            MinimizerType::ReducedMinimizer { k0 } => {
                density(text, l, |lmer| reduced_minimizer(lmer, k, *k0))
            }
        }
    }

    /// For BdAnchor and ReducedMinimizer, try r/k0 from 3 to 7.
    /// For Miniception, try k0 from 3 to 7, and k0 = k-w = k-(l-k+1) = 2k-l-1.
    fn try_params(&self, l: usize, k: usize) -> Vec<Self> {
        match self {
            MinimizerType::Minimizer
            | MinimizerType::RobustMinimizer
            | MinimizerType::BiMinimizer
            | MinimizerType::BiMinimizerBot => vec![*self],
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
            MinimizerType::ReducedMinimizer { .. } => {
                // k <= (l+k0+1)/2
                // 2k <= l + k0 + 1
                // 2k - l - 1 <= k0
                let k0_min = 1.max((2 * k - 1).saturating_sub(l));
                let k0_max = 10.max(k0_min + 2).min(k);
                (k0_min..=k0_max)
                    .map(|k0| MinimizerType::ReducedMinimizer { k0 })
                    .collect()
            }
        }
    }
}

#[derive(clap::Subcommand)]
enum Command {
    Run {
        /// Test kmers up to l.
        #[clap(short, default_value_t = 40)]
        l: usize,
        /// Minimizer length.
        #[arg(short, default_value_t = 10)]
        k: usize,
        #[clap(subcommand)]
        tp: MinimizerType,
    },
    Eval,
}

/// Print expected and measured densities for different l and r.
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
        Command::Run { tp, l, k } => {
            eprintln!("Running {tp:?}:");
            let d = tp.density(text, l, k);
            eprintln!("  Density: {:.3}", d);
        }
        Command::Eval => {
            let base_types = [
                MinimizerType::Minimizer,
                MinimizerType::RobustMinimizer,
                MinimizerType::BdAnchor { r: 0 },
                MinimizerType::Miniception { k0: 0 },
                MinimizerType::BiMinimizer,
                MinimizerType::BiMinimizerBot,
                MinimizerType::ReducedMinimizer { k0: 0 },
            ];

            let mut results = vec![];

            for k in [
                4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128,
            ] {
                for w in [8, 32, 128] {
                    let l = k + w - 1;
                    for tp in base_types.iter() {
                        let tps = &tp.try_params(l, k);
                        let (d, tp) = tps
                            .iter()
                            .map(|tp| (tp.density(text, l, k), tp))
                            .min_by(|&(ld, _), &(rd, _)| {
                                if ld < rd {
                                    std::cmp::Ordering::Less
                                } else {
                                    std::cmp::Ordering::Greater
                                }
                            })
                            .unwrap();
                        results.push(Result {
                            sigma: args.sigma,
                            k,
                            w,
                            l,
                            tp: *tp,
                            density: d,
                        });
                        eprintln!("k={k} w={w} l={l} tp={tp:?} d={d:.3}",);
                    }
                }
            }

            let result_json = serde_json::to_string(&results).unwrap();
            println!("{}", result_json);
        }
    }
}
