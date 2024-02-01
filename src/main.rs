#![feature(exclusive_range_pattern, type_alias_impl_trait)]
use std::{
    cmp::Reverse, collections::VecDeque, f32::consts::PI, io::Write, iter::zip, path::PathBuf,
    sync::Mutex,
};

use clap::Parser;
use itertools::Itertools;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
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

fn minimizer(s: &[u8], k: usize) -> usize {
    assert!(k > 0);
    assert!(k <= s.len());
    s.windows(k)
        .enumerate()
        .min_by_key(|&(i, kmer)| (h(kmer), Reverse(i)))
        .unwrap()
        .0
}

fn text_minimizers(text: &[u8], l: usize, k: usize) -> Vec<usize> {
    let mut q = IQ::new();
    let w = l - k + 1;
    let mut kmers = text.windows(k).enumerate();
    for (j, kmer) in kmers.by_ref().take(w - 1) {
        q.push(j, h(kmer));
    }
    // i: position of lmer
    // j: position of kmer
    kmers
        .enumerate()
        .map(|(i, (j, kmer))| {
            q.push(j, h(kmer));
            q.pop(i).unwrap().0
        })
        .dedup()
        .collect()
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

fn text_miniception(text: &[u8], l: usize, k: usize, k0: usize) -> Vec<usize> {
    let w = l - k + 1;
    // The number of k0-mers in a kmer.
    let w0 = k - k0;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    // Queue of all k0-mers.
    let mut q0 = IQ::new();
    // Queue of filtered k-mers.
    let mut q = IQ::new();

    // i: position of lmer
    // j: position of kmer
    // j0: position of k0mer

    // 1: init k0-mers.
    let mut k0mers = text.windows(k0).enumerate();
    for (j0, k0mer) in k0mers.by_ref().take(w0) {
        q0.push(j0, h(k0mer));
    }

    // 2: init k-mers.
    let mut kmers = text.windows(k).enumerate().zip(k0mers);
    for ((j, kmer), (j0, k0mer)) in kmers.by_ref().take(w - 1) {
        q0.push(j0, h(k0mer));
        let min_pos = q0.pop(j).unwrap().0;
        if min_pos == j || min_pos == j + w0 {
            q.push(j, h(kmer));
        }
    }

    // 3: Iterate l-mers.
    kmers
        .enumerate()
        .map(|(i, ((j, kmer), (j0, k0mer)))| {
            q0.push(j0, h(k0mer));
            let min_pos = q0.pop(j).unwrap().0;
            if min_pos == j || min_pos == j + w0 {
                q.push(j, h(kmer));
            }

            q.pop(i).unwrap().0
        })
        .dedup()
        .collect()
}

/// Sort filtered kmers by:
/// - first prefixes by h(k0), then suffixes by h(k0)
///
/// For small k, other orderings are worse, such as:
/// - first prefixes by h(k0), then suffixes by -h(k0)
/// - first prefixes by -h(k0), then suffixes by h(k0)
/// - first prefixes by -h(k0), then suffixes by -h(k0)
///
/// Sorting by h(k0) directly is same as mod_minimizers and best for large k.
fn miniception_new(s: &[u8], k: usize, k0: usize) -> usize {
    let w = s.len() - k + 1;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    let w0 = k - k0;
    let _l0 = w0 + k0 - 1;

    s.windows(k)
        .enumerate()
        .filter_map(|(i, kmer)| {
            let j = minimizer(kmer, k0);
            assert!(j <= k - k0);
            if j == 0 || j == k - k0 {
                Some((i, kmer, j == 0, h(&kmer[j..j + k0])))
            } else {
                None
            }
        })
        .min_by_key(|&(i, w, b, hk0)| (b, hk0, h(w), Reverse(i)))
        .unwrap()
        .0
}

// TODO: dedup
fn text_miniception_new(text: &[u8], l: usize, k: usize, k0: usize) -> Vec<usize> {
    let w = l - k + 1;
    // The number of k0-mers in a kmer.
    let w0 = k - k0;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    // Queue of all k0-mers.
    let mut q0 = IQ::new();
    // Queue of filtered k-mers.
    let mut q = IQ::new();

    // i: position of lmer
    // j: position of kmer
    // j0: position of k0mer

    // 1: init k0-mers.
    let mut k0mers = text.windows(k0).enumerate();
    for (j0, k0mer) in k0mers.by_ref().take(w0) {
        q0.push(j0, h(k0mer));
    }

    // 2: init k-mers.
    let mut kmers = text.windows(k).enumerate().zip(k0mers);
    for ((j, kmer), (j0, k0mer)) in kmers.by_ref().take(w - 1) {
        q0.push(j0, h(k0mer));
        let min_pos = q0.pop(j).unwrap().0;
        assert!(j <= min_pos && min_pos <= j + w0);
        if min_pos == j || min_pos == j + w0 {
            // TODO: Test without h(kmer)?
            q.push(j, (min_pos == j, h(&text[min_pos..min_pos + k0]), h(kmer)));
        }
    }

    // 3: Iterate l-mers.
    kmers
        .enumerate()
        .map(|(i, ((j, kmer), (j0, k0mer)))| {
            q0.push(j0, h(k0mer));
            let min_pos = q0.pop(j).unwrap().0;
            if min_pos == j || min_pos == j + w0 {
                q.push(j, (min_pos == j, h(&text[min_pos..min_pos + k0]), h(kmer)));
            }

            q.pop(i).unwrap().0
        })
        .dedup()
        .collect()
}

fn robust_biminimizer(s: &[u8], k: usize, last: &mut usize) -> usize {
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

/// Find minimal t-mer at pos idx. Then select idx % w.
fn lr_minimizer(s: &[u8], k: usize, t: usize) -> usize {
    let l = s.len();
    let w = l - k + 1;
    let idx = s
        .windows(t)
        .enumerate()
        .min_by_key(|&(i, w)| (h(w), Reverse(i)))
        .unwrap()
        .0;
    let i = if idx >= w { idx - w } else { idx };
    assert!(i + k <= l);
    i
}

fn text_lr_minimizers(text: &[u8], l: usize, k: usize, t: usize) -> Vec<usize> {
    let mut q = IQ::new();
    let w = l - k + 1;
    let wt = l - t + 1;
    let mut tmers = text.windows(t).enumerate();

    for (j, tmer) in tmers.by_ref().take(wt - 1) {
        q.push(j, h(tmer));
    }
    // i: position of lmer
    // j: position of tmer
    tmers
        .enumerate()
        .map(|(i, (j, tmer))| {
            q.push(j, h(tmer));
            let idx = q.pop(i).unwrap().0 - i;
            i + if idx >= w { idx - w } else { idx }
        })
        .dedup()
        .collect()
}

/// Find minimal t-mer at pos idx. Then select idx % w.
fn mod_minimizer(s: &[u8], k: usize, t: usize) -> usize {
    let l = s.len();
    let w = l - k + 1;
    let idx = s
        .windows(t)
        .enumerate()
        .min_by_key(|&(i, w)| (h(w), Reverse(i)))
        .unwrap()
        .0;
    idx % w
}

/// FastMod32, using the low 32 bits of the hash.
/// Taken from https://github.com/lemire/fastmod/blob/master/include/fastmod.h
#[derive(Copy, Clone, Debug)]
pub struct FM32 {
    d: u64,
    m: u64,
}
impl FM32 {
    fn new(d: usize) -> Self {
        assert!(d <= u32::MAX as usize);
        Self {
            d: d as u64,
            m: u64::MAX / d as u64 + 1,
        }
    }
    fn reduce(self, h: usize) -> usize {
        let lowbits = self.m * (h as u64);
        ((lowbits as u128 * self.d as u128) >> 64) as usize
    }
}

fn text_mod_minimizers(text: &[u8], l: usize, k: usize, t: usize) -> Vec<usize> {
    let mut q = IQ::new();
    let w = l - k + 1;
    let wt = l - t + 1;
    let mut tmers = text.windows(t).enumerate();

    let fastmod_w = FM32::new(w);

    for (j, tmer) in tmers.by_ref().take(wt - 1) {
        q.push(j, h(tmer));
    }
    // i: position of lmer
    // j: position of tmer
    tmers
        .enumerate()
        .map(|(i, (j, tmer))| {
            q.push(j, h(tmer));
            i + fastmod_w.reduce(q.pop(i).unwrap().0 - i)
        })
        .dedup()
        .collect()
}

/// Asymptotic (in k) optimal minimizers:
/// - Assume k=x*w; sum all i mod w positions.
/// - Take minimizer if first coordinate is max or at most sigma away from max.
/// - If multiple, take random.
///
/// - Ignore for k < w.
/// - Ignore trailing chars.
///
/// Simplified:
/// - For all k-mers compute the rotational sum of k/w w-strided values.
/// - Take k-mer that maximizes the first coordinate.
fn rot_minimizer(s: &[u8], k: usize) -> usize {
    let w = s.len() - k + 1;
    assert!(k >= w);
    (0..=s.len() - k)
        .max_by_key(|&i| {
            s[i..i + k]
                .iter()
                .step_by(w)
                .map(|&x| x as u64)
                .sum::<u64>()
        })
        .unwrap()
}

fn decycling_minimizer_init(k: usize) -> Vec<f32> {
    // (0..k).map(|i| (i as f32).sin()).collect()
    (0..k)
        .map(|i| (i as f32 / k as f32 * 2. * PI).sin())
        .collect()
}

/// Check if this is in the decyling set.
/// Note that we ignore the edge case where ix=0 and ix2=0.
fn is_in_decycling_set(kmer: &[u8], cs: &Vec<f32>) -> bool {
    let ix = zip(kmer, cs).map(|(&x, c)| x as f32 * c).sum::<f32>();
    let ix2 = *kmer.last().unwrap() as f32 * cs[0]
        + zip(kmer, &cs[1..]).map(|(&x, c)| x as f32 * c).sum::<f32>();
    ix > 0. && ix2 <= 0.
}
fn is_in_double_decycling_set(kmer: &[u8], cs: &Vec<f32>) -> bool {
    let ix = zip(kmer, cs).map(|(&x, c)| x as f32 * c).sum::<f32>();
    let ix2 = *kmer.last().unwrap() as f32 * cs[0]
        + zip(kmer, &cs[1..]).map(|(&x, c)| x as f32 * c).sum::<f32>();
    ix < 0. && ix2 >= 0.
}

/// Decycling minimizer using Mykkelvelt embedding
fn decycling_minimizer(s: &[u8], k: usize, cs: &Vec<f32>) -> usize {
    s.windows(k)
        .enumerate()
        .min_by_key(|&(i, kmer)| (!is_in_decycling_set(kmer, cs), h(kmer), Reverse(i)))
        .unwrap()
        .0
}

/// Decycling minimizer using Mykkelvelt embedding
fn double_decycling_minimizer(s: &[u8], k: usize, cs: &Vec<f32>) -> usize {
    s.windows(k)
        .enumerate()
        .min_by_key(|&(i, kmer)| {
            (
                !is_in_decycling_set(kmer, cs),
                !is_in_double_decycling_set(kmer, cs),
                h(kmer),
                Reverse(i),
            )
        })
        .unwrap()
        .0
}

fn collect_anchors(text: &[u8], l: usize, mut scheme: impl FnMut(&[u8]) -> usize) -> Vec<usize> {
    text.windows(l)
        .enumerate()
        .map(|(i, w)| {
            let x = i + scheme(w);
            // eprintln!(" {x:>4}");
            x
        })
        .dedup()
        .collect_vec()
}

/// Compute the density of the sampling scheme.
/// A function to select the minimal kmer from an lmer.
/// The length of the slice is l.
/// k must be bound by the function.
///
/// Must return a value in [0, l - k].
fn density(text: &[u8], l: usize, scheme: impl FnMut(&[u8]) -> usize) -> f64 {
    let mut anchors = collect_anchors(text, l, scheme);
    // Dedup anchors to ensure valid results for non-forward schemes (bd-anchors).
    // TODO: Analyze non-forward schemes.
    anchors.sort();
    anchors.dedup();

    anchors.len() as f64 / text.windows(l).len() as f64
}

#[derive(Clone, Copy, clap::Subcommand, Debug, Serialize)]
#[serde(tag = "minimizer_type")]
enum MinimizerType {
    Minimizer,
    BdAnchor { r: usize },
    Miniception { k0: usize },
    MiniceptionNew { k0: usize },
    BiMinimizer,
    ModMinimizer { k0: usize },
    LrMinimizer { k0: usize },
    RotMinimizer,
    DecyclingMinimizer,
    DoubleDecyclingMinimizer,
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
    #[inline(never)]
    fn density(&self, text: &[u8], l: usize, k: usize) -> f64 {
        let d = |minis: Vec<usize>| minis.len() as f64 / text.windows(l).len() as f64;
        match self {
            MinimizerType::Minimizer => d(text_minimizers(text, l, k)),
            MinimizerType::BdAnchor { r } => density(text, l, |lmer| bd_anchor(lmer, *r)),
            MinimizerType::Miniception { k0 } => d(text_miniception(text, l, k, *k0)),
            MinimizerType::MiniceptionNew { k0 } => d(text_miniception_new(text, l, k, *k0)),
            MinimizerType::BiMinimizer => {
                let last = &mut 0;
                density(text, l, move |lmer| robust_biminimizer(lmer, k, last))
            }
            MinimizerType::ModMinimizer { k0 } => d(text_mod_minimizers(text, l, k, *k0)),
            MinimizerType::LrMinimizer { k0 } => d(text_lr_minimizers(text, l, k, *k0)),
            MinimizerType::RotMinimizer => density(text, l, move |lmer| rot_minimizer(lmer, k)),
            MinimizerType::DecyclingMinimizer => {
                let cs = decycling_minimizer_init(k);
                density(text, l, move |lmer| decycling_minimizer(lmer, k, &cs))
            }
            MinimizerType::DoubleDecyclingMinimizer => {
                let cs = decycling_minimizer_init(k);
                density(text, l, move |lmer| {
                    double_decycling_minimizer(lmer, k, &cs)
                })
            }
        }
    }

    fn try_params(&self, l: usize, k: usize) -> Vec<Self> {
        let w = l - k + 1;
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
        /// Test kmers up to l.
        #[clap(short, default_value_t = 40)]
        l: usize,
        /// Minimizer length.
        #[arg(short, default_value_t = 10)]
        k: usize,
        #[clap(subcommand)]
        tp: MinimizerType,
    },
    Eval {
        /// Write to file.
        #[clap(short, long)]
        output: Option<PathBuf>,
    },
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
        Command::Eval { output } => {
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

            let results = Mutex::new(vec![]);

            let ks = [
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 52, 56, 60, 64,
            ];
            let ws = [8, 16, 32];
            let kws = ks.into_iter().cartesian_product(ws).collect_vec();
            kws.par_iter().for_each(|&(k, w)| {
                let l = k + w - 1;
                for tp in base_types.iter() {
                    let tps = &tp.try_params(l, k);
                    if tps.is_empty() {
                        continue;
                    }
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
                    results.lock().unwrap().push(Result {
                        sigma: args.sigma,
                        k,
                        w,
                        l,
                        tp: *tp,
                        density: d,
                    });
                    eprintln!("k={k} w={w} l={l} tp={tp:?} d={d:.3}",);
                }
            });

            if let Some(output) = output {
                let result_json = serde_json::to_string(&results).unwrap();
                let mut file = std::fs::File::create(output).unwrap();
                file.write_all(result_json.as_bytes()).unwrap();
            }
        }
    }
}

/// For increasing timestamps, store increasing values.
struct IQ<K: Ord> {
    q: VecDeque<(usize, K)>,
}

impl<K: Ord> IQ<K>
where
    K: Copy,
{
    fn new() -> Self {
        Self { q: VecDeque::new() }
    }

    fn push(&mut self, t: usize, k: K) {
        while let Some(&back) = self.q.back() {
            if back.1 >= k {
                self.q.pop_back();
            } else {
                break;
            }
        }
        self.q.push_back((t, k));
    }

    /// Return the minimal element with timestamp >= t.
    fn pop(&mut self, t: usize) -> Option<&(usize, K)> {
        while let Some(&front) = self.q.front() {
            if front.0 < t {
                self.q.pop_front();
            } else {
                break;
            }
        }
        self.q.front()
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
                let l = k + w - 1;
                let anchors = collect_anchors(&text, l, |lmer| minimizer(lmer, k));
                let minimizers = text_minimizers(&text, l, k);
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
                    let anchors = collect_anchors(&text, l, |lmer| mod_minimizer(lmer, k, t));
                    let minimizers = text_mod_minimizers(&text, l, k, t);
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
                    let anchors = collect_anchors(&text, l, |lmer| lr_minimizer(lmer, k, t));
                    let minimizers = text_lr_minimizers(&text, l, k, t);
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
                let l = k + w - 1;
                for k0 in k.saturating_sub(w).max(1)..=k {
                    let anchors = collect_anchors(&text, l, |lmer| super::miniception(lmer, k, k0));
                    let minimizers = text_miniception(&text, l, k, k0);
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
                let l = k + w - 1;
                for k0 in k.saturating_sub(w).max(1)..=k {
                    let anchors =
                        collect_anchors(&text, l, |lmer| super::miniception_new(lmer, k, k0));
                    let minimizers = text_miniception_new(&text, l, k, k0);
                    assert_eq!(anchors, minimizers);
                }
            }
        }
    }
}
