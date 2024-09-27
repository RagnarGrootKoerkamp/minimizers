// TODO: Implement fast version of AntiLex sort, skipping NtHash.
// TODO: Make the 'linearized' `minimizers_simd_it` and `nthash32_simd_it` iterators faster.
//       Currently they are 4x slower than the `*_par_it` versions.

pub mod bruteforce;
pub mod cycle_partition_lp;
pub mod de_bruijn_seq;
pub mod local_scheme;
pub mod monotone_queue;
pub mod order;
#[cfg(feature = "python-bindings")]
mod py;
pub mod simd;

use crate::monotone_queue::MonotoneQueue;
use itertools::Itertools;
use num::{complex::Complex64, Zero};
use order::*;
use rand_chacha::{
    rand_core::{RngCore, SeedableRng},
    ChaChaRng,
};
use serde::{Deserialize, Serialize};
use std::{
    cell::Cell,
    cmp::{max, Reverse},
    f64::consts::PI,
};

/// An iterator over *all* minimizer positions. Not yet deduplicated.
///
/// NOTE: For non-forward schemes, positions may be returned twice.
pub trait MinimizerIt: Iterator<Item = usize> {}
impl<I: Iterator<Item = usize>> MinimizerIt for I {}

pub trait SamplingScheme {
    fn w(&self) -> usize {
        panic!();
    }
    fn k(&self) -> usize {
        panic!();
    }
    fn l(&self) -> usize;

    /// Sample a single lmer.
    fn sample(&self, _lmer: &[u8]) -> usize {
        unimplemented!("Use stream() instead.");
    }

    /// Sample all lmers in a text.
    /// This default implementation simply calls `sample` on each lmer.
    // TODO: Take an iterator over u8 instead?
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        self.stream_naive(text)
    }

    /// Sample all lmers in a cyclic text of length `len`.
    /// Text must have length at least `len + l-1`, and the additional
    /// characters must equal a prefix of the text.
    ///
    /// Returns the number of selected minimizers.
    fn cyclic_text_density(&self, text: &[u8], len: usize) -> usize {
        let text = &text[..len + self.l() - 1];
        let mut poss: Vec<_> = self.stream(text).iter().map(|p| p % len).collect();
        poss.sort();
        poss.dedup();
        poss.len()
    }

    fn is_forward(&self, text: &[u8]) -> bool {
        self.stream(text)
            .iter()
            .tuple_windows()
            .all(|(a, b)| a <= b)
    }

    /// Sample all lmers in a text.
    /// This default implementation simply calls `sample` on each lmer.
    fn stream_naive(&self, text: &[u8]) -> Vec<usize> {
        text.windows(self.l())
            .enumerate()
            .map(|(i, lmer)| i + self.sample(lmer))
            .collect()
    }
}
/// Generate a random string.
pub fn generate_random_string(n: usize, sigma: usize) -> Vec<u8> {
    let mut rng = ChaChaRng::seed_from_u64(213456);
    (0..n)
        .map(|_| (((rng.next_u64() as usize) % sigma) as u8))
        .collect()
}

/// Returns:
/// - density
/// - position distribution
/// - distance distribution
/// - transfer distribution
pub fn collect_stats(
    w: usize,
    text: &[u8],
    scheme: &dyn SamplingScheme,
) -> (f64, Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let it = scheme.stream(text);

    let mut n = 0;
    let mut anchors = 0;
    let mut ps = vec![0; w];
    let mut ds = vec![0; 2 * w + 1];
    let mut transfer = vec![vec![0; w]; w];
    let mut last = 0;
    for (i, idx) in it.into_iter().enumerate() {
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

/// Compute statistics on number of sampled positions on cycles of a given length.
pub fn cycle_stats(l: usize, text: &[u8], scheme: &dyn SamplingScheme) -> (f64, Vec<f64>) {
    let mut cycle = vec![0; 4 * l];
    let mut stats = vec![0; l + 1];
    let mut total = 0;
    let mut num_windows = 0;
    let mut num_cycles = 0;

    // Find length-l cycles.
    for c in text.chunks_exact(l) {
        // fill cycle with c twice
        cycle[..l].copy_from_slice(c);
        cycle[l..2 * l].copy_from_slice(c);
        cycle[2 * l..3 * l].copy_from_slice(c);
        cycle[3 * l..4 * l].copy_from_slice(c);
        let samples = scheme.cyclic_text_density(&cycle, l);
        stats[samples] += 1;
        total += samples;
        num_windows += l;
        num_cycles += 1;
    }
    let density = total as f64 / num_windows as f64;
    let stats = stats
        .into_iter()
        .map(|c| c as f64 / num_cycles as f64)
        .collect();
    (density, stats)
}

#[derive(Clone, Copy, clap::Subcommand, Debug, Serialize, Deserialize)]
#[serde(tag = "minimizer_type")]
pub enum MinimizerType {
    Minimizer {
        ao: bool,
    },
    BdAnchor {
        r: usize,
    },
    SusAnchor {
        ao: bool,
        modulo: bool,
    },
    Miniception {
        k0: usize,
        ao: bool,
        aot: bool,
    },
    MiniceptionNew {
        k0: usize,
    },
    ModSampling {
        k0: usize,
    },
    LrMinimizer,
    ModMinimizer {
        r: usize,
        aot: bool,
    },
    RotMinimizer,
    AltRotMinimizer,
    DecyclingMinimizer,
    DoubleDecyclingMinimizer {
        ao: bool,
    },
    Bruteforce,
    OpenSyncmerMinimizer {
        t: usize,
    },
    ClosedSyncmerMinimizer {
        t: usize,
        h: usize,

        loose: bool,
        open: bool,
    },
    OpenClosedSyncmerMinimizer {
        t: usize,
    },
    FracMin {
        f: usize,
    },
    OcModMinimizer {
        t: usize,
        offset: usize,
        use_closed: bool,
        prefer_prefix: bool,
        open_tmer: bool,
        closed_tmer: bool,
        other_tmer: bool,
        ao: bool,
        aot: bool,
    },
}

/// TODO: Analyze non-forward schemes.
impl MinimizerType {
    pub fn build(&self, w: usize, k: usize, sigma: usize) -> Box<dyn SamplingScheme> {
        let o = RandomOrder;
        let al = AntiLex;
        match self {
            MinimizerType::Minimizer { ao } => {
                if !ao {
                    Box::new(Minimizer::new(k, w, o))
                } else {
                    Box::new(Minimizer::new(k, w, al))
                }
            }
            MinimizerType::BdAnchor { r } => {
                assert_eq!(k, 1);
                Box::new(BdAnchor::new(w, *r))
            }
            MinimizerType::SusAnchor { ao, modulo } => {
                if !ao {
                    Box::new(SusAnchor::new(w, k, Lex, *modulo))
                } else {
                    Box::new(SusAnchor::new(w, k, al, *modulo))
                }
            }
            MinimizerType::Miniception { k0, ao, aot } => {
                if !ao {
                    if !aot {
                        Box::new(Miniception::new(w, k, *k0, o, o))
                    } else {
                        Box::new(Miniception::new(w, k, *k0, o, al))
                    }
                } else {
                    if !aot {
                        Box::new(Miniception::new(w, k, *k0, al, o))
                    } else {
                        Box::new(Miniception::new(w, k, *k0, al, al))
                    }
                }
            }
            MinimizerType::MiniceptionNew { k0 } => Box::new(MiniceptionNew::new(w, k, *k0, o)),
            MinimizerType::ModSampling { k0 } => Box::new(ModSampling::new(k, w, *k0, o)),
            MinimizerType::LrMinimizer => Box::new(ModSampling::lr_minimizer(k, w, o)),
            MinimizerType::ModMinimizer { r, aot } => {
                if !aot {
                    Box::new(ModSampling::mod_minimizer(k, w, *r, o))
                } else {
                    Box::new(ModSampling::mod_minimizer(k, w, *r, al))
                }
            }
            MinimizerType::RotMinimizer => Box::new(RotMinimizer::new(k, w, o)),
            MinimizerType::AltRotMinimizer => Box::new(AltRotMinimizer::new(k, w, o)),
            MinimizerType::DecyclingMinimizer => Box::new(Decycling::new(k, w, o, false)),
            MinimizerType::DoubleDecyclingMinimizer { ao } => {
                if !ao {
                    Box::new(Decycling::new(k, w, o, true))
                } else {
                    Box::new(Decycling::new(k, w, al, true))
                }
            }
            MinimizerType::Bruteforce => {
                let m = bruteforce::bruteforce_minimizer(k, w, sigma).1;
                Box::new(m)
            }
            MinimizerType::OpenSyncmerMinimizer { t } => {
                Box::new(OpenSyncmer::new(k, w, *t, true, false))
            }
            MinimizerType::ClosedSyncmerMinimizer { t, loose, open, h } => {
                Box::new(ClosedSyncmer::new(k, w, *t, *loose, *open, *h))
            }
            MinimizerType::OpenClosedSyncmerMinimizer { t } => {
                Box::new(OpenSyncmer::new(k, w, *t, true, true))
            }
            MinimizerType::FracMin { f } => Box::new(FracMin::new(k, w, *f)),
            MinimizerType::OcModMinimizer {
                t,
                offset,
                use_closed,
                prefer_prefix,
                open_tmer,
                closed_tmer,
                other_tmer,
                ao,
                aot,
            } => {
                if !ao {
                    if !aot {
                        Box::new(OcModMinimizer::new(
                            k,
                            w,
                            *t,
                            *offset,
                            *use_closed,
                            *prefer_prefix,
                            *open_tmer,
                            *closed_tmer,
                            *other_tmer,
                            o,
                            o,
                        ))
                    } else {
                        Box::new(OcModMinimizer::new(
                            k,
                            w,
                            *t,
                            *offset,
                            *use_closed,
                            *prefer_prefix,
                            *open_tmer,
                            *closed_tmer,
                            *other_tmer,
                            o,
                            al,
                        ))
                    }
                } else {
                    if !aot {
                        Box::new(OcModMinimizer::new(
                            k,
                            w,
                            *t,
                            *offset,
                            *use_closed,
                            *prefer_prefix,
                            *open_tmer,
                            *closed_tmer,
                            *other_tmer,
                            al,
                            o,
                        ))
                    } else {
                        Box::new(OcModMinimizer::new(
                            k,
                            w,
                            *t,
                            *offset,
                            *use_closed,
                            *prefer_prefix,
                            *open_tmer,
                            *closed_tmer,
                            *other_tmer,
                            al,
                            al,
                        ))
                    }
                }
            }
        }
    }

    pub fn try_params(&self, w: usize, k: usize) -> Vec<Self> {
        let l = w + k - 1;
        match self {
            MinimizerType::BdAnchor { .. } => {
                let r_max = k;
                (0.min(r_max)..=10.min(r_max))
                    .map(|r| MinimizerType::BdAnchor { r })
                    .collect()
            }
            MinimizerType::Miniception { ao, aot, .. } => {
                let k0 = max(k as isize - w as isize, 4) as usize;
                return vec![MinimizerType::Miniception {
                    k0,
                    ao: *ao,
                    aot: *aot,
                }];
            }
            MinimizerType::MiniceptionNew { .. } => {
                let k0 = max(k as isize - w as isize, 4) as usize;
                return vec![MinimizerType::MiniceptionNew { k0 }];
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
                vec![MinimizerType::OpenSyncmerMinimizer { t: 4 }]
            }
            MinimizerType::OpenClosedSyncmerMinimizer { .. } => {
                vec![MinimizerType::OpenClosedSyncmerMinimizer { t: 4 }]
            }
            MinimizerType::FracMin { .. } => (1..w).map(|f| MinimizerType::FracMin { f }).collect(),
            _ => {
                vec![*self]
            }
        }
    }
}

/// Classic minimizers with respect to some order on kmers `O`.
pub struct Minimizer<O: DirectedOrder> {
    k: usize,
    w: usize,
    l: usize,
    o: O,
}

impl<O: DirectedOrder> Minimizer<O> {
    pub fn new(k: usize, w: usize, o: O) -> Self {
        assert!(k > 0);
        let l = k + w - 1;
        Self { k, w, l, o }
    }
    pub fn ord(&self) -> &O {
        &self.o
    }
}

impl<O: DirectedOrder> SamplingScheme for Minimizer<O> {
    fn w(&self) -> usize {
        self.w
    }
    fn k(&self) -> usize {
        self.k
    }
    fn l(&self) -> usize {
        self.l
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        debug_assert_eq!(lmer.len(), self.l);
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|&(_i, kmer)| self.o.key(kmer))
            .unwrap()
            .0
    }

    // TODO: Rolling hash using NtHash.
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        let mut q = monotone_queue::MonotoneQueue::new();
        let mut kmers = text.windows(self.k).enumerate();
        // Push the first w-1 kmers onto the queue.
        for (j, kmer) in kmers.by_ref().take(self.w - 1) {
            q.push(j, self.o.key(kmer));
        }
        // i: position of lmer
        // j: position of kmer
        kmers
            .enumerate()
            .map(move |(i, (j, kmer))| {
                q.push(j, self.o.key(kmer));
                q.pop(i).unwrap().0
            })
            .collect()
    }
}

/// Classic random minimizers.
pub type RandomMinimizer = Minimizer<RandomOrder>;

#[derive(Debug, Clone)]
pub struct ExplicitLocalScheme {
    pub k: usize,
    pub w: usize,
    pub sigma: usize,
    /// The index in [w] to choose for each of sigma^(k+w-1) possible l-mers.
    pub map: Vec<u8>,
}

impl SamplingScheme for ExplicitLocalScheme {
    fn w(&self) -> usize {
        self.w
    }
    fn k(&self) -> usize {
        self.k
    }
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        debug_assert_eq!(lmer.len(), self.l());
        debug_assert_eq!(self.map.len(), self.sigma.pow(self.l() as u32));
        let mut v = 0;
        for c in lmer {
            assert!(*c < self.sigma as u8);
            v = self.sigma * v + *c as usize;
        }
        self.map[v] as usize
    }
}

pub struct BdAnchor {
    w: usize,
    r: usize,
}

impl BdAnchor {
    pub fn new(w: usize, r: usize) -> Self {
        Self { w, r }
    }
}

impl SamplingScheme for BdAnchor {
    fn l(&self) -> usize {
        self.w
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let w = self.w;
        debug_assert_eq!(lmer.len(), w);
        let mut best = 0;
        for i in 1..w.saturating_sub(self.r) {
            for j in 0..w {
                if lmer[(i + j) % w] != lmer[(best + j) % w] {
                    if lmer[(i + j) % w] < lmer[(best + j) % w] {
                        best = i;
                    }
                    break;
                }
            }
        }
        best
    }
}

/// NOTE: O should be Lex or AntiLex order. Random order will not be good.
pub struct SusAnchor<O: Order> {
    w: usize,
    k: usize,
    o: O,
    modulo: bool,
}

impl<O: Order> SusAnchor<O> {
    pub fn new(w: usize, k: usize, o: O, modulo: bool) -> Self {
        Self { w, k, o, modulo }
    }
}

impl<O: Order> SamplingScheme for SusAnchor<O> {
    fn l(&self) -> usize {
        self.w + self.k - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let mut best = (self.o.key(lmer), 0);
        if self.modulo {
            if lmer.iter().all(|&c| c == 0) {
                return 1;
            }
            let t = (self.k - 1) % self.w + 1;
            for i in 1..=self.l() - t {
                best = best.min((self.o.key(&lmer[i..]), i));
            }
            best.1 % self.w
        } else {
            for i in 1..self.w {
                best = best.min((self.o.key(&lmer[i..]), i));
            }
            best.1
        }
    }
}

pub struct Miniception<O: Order, OT: DirectedOrder> {
    w: usize,
    k: usize,
    l: usize,
    k0: usize,
    w0: usize,
    /// The order used for k-mers.
    o: O,
    /// The order used for k0-mers.
    ot: Minimizer<OT>,
}

impl<O: Order, OT: DirectedOrder> Miniception<O, OT> {
    pub fn new(w: usize, k: usize, k0: usize, o: O, ot: OT) -> Self {
        assert!(k0 >= k.saturating_sub(w));
        assert!(k0 <= k);
        let w0 = k - k0;
        Self {
            w,
            k,
            l: k + w - 1,
            k0,
            w0,
            o,
            ot: Minimizer::new(k0, k - k0 + 1, ot),
        }
    }
}

impl<O: Order, OT: Order> SamplingScheme for Miniception<O, OT> {
    fn l(&self) -> usize {
        self.l
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .filter(|(_, kmer)| {
                let j = self.ot.sample(kmer);
                j == 0 || j == self.k - self.k0
            })
            .min_by_key(|&(_i, kmer)| self.o.key(kmer))
            .unwrap()
            .0
    }

    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        // Queue of all k0-mers.
        let mut q0 = MonotoneQueue::new();
        // Queue of filtered k-mers.
        let mut q = MonotoneQueue::new();

        // i: position of lmer
        // j: position of kmer
        // j0: position of k0mer

        let o0 = &self.ot.o;

        // 1: init k0-mers.
        let mut k0mers = text.windows(self.k0).enumerate();
        for (j0, k0mer) in k0mers.by_ref().take(self.w0) {
            q0.push(j0, Order::key(o0, k0mer));
        }

        // 2: init k-mers.
        let mut kmers = text.windows(self.k).enumerate().zip(k0mers);
        for ((j, kmer), (j0, k0mer)) in kmers.by_ref().take(self.w - 1) {
            q0.push(j0, Order::key(o0, k0mer));
            let min_pos = q0.pop(j).unwrap().0;
            if min_pos == j || min_pos == j + self.w0 {
                q.push(j, self.o.key(kmer));
            }
        }

        // 3: Iterate l-mers.
        kmers
            .enumerate()
            .map(move |(i, ((j, kmer), (j0, k0mer)))| {
                q0.push(j0, Order::key(o0, k0mer));
                let min_pos = q0.pop(j).unwrap().0;
                if min_pos == j || min_pos == j + self.w0 {
                    q.push(j, self.o.key(kmer));
                }

                q.pop(i).unwrap().0
            })
            .collect()
    }
}

/// Modified version of miniception with a different sort order on kmer.
/// Instead of sorting them by h(s), the are ordered by:
/// - First those for which the minimal k0-mer is a prefix, sorted by h(t).
/// - Then those for which the minimal k0-mer is a suffix, sorted by h(t).
pub struct MiniceptionNew<O: Order> {
    w: usize,
    k: usize,
    l: usize,
    k0: usize,
    w0: usize,
    /// The order used for k-mers.
    o: O,
    /// The order used for k0-mers.
    rand_mini: RandomMinimizer,
}

impl<O: Order> MiniceptionNew<O> {
    pub fn new(w: usize, k: usize, k0: usize, o: O) -> Self {
        assert!(k0 >= k.saturating_sub(w));
        assert!(k0 <= k);
        let w0 = k - k0;
        Self {
            w,
            k,
            l: k + w - 1,
            k0,
            w0,
            o,
            rand_mini: RandomMinimizer::new(k0, k - k0 + 1, RandomOrder),
        }
    }
}

impl<O: Order> SamplingScheme for MiniceptionNew<O> {
    fn l(&self) -> usize {
        self.l
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let o0 = &self.rand_mini.o;
        lmer.windows(self.k)
            .enumerate()
            .filter_map(|(i, kmer)| {
                let j = self.rand_mini.sample(kmer);
                if j == 0 || j == self.k - self.k0 {
                    Some((i, kmer, j == 0, Order::key(o0, &kmer[j..j + self.k0])))
                } else {
                    None
                }
            })
            .min_by_key(|&(_i, w, b, hk0)| (b, hk0, self.o.key(w)))
            .unwrap()
            .0
    }

    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        // Queue of all k0-mers.
        let mut q0 = MonotoneQueue::new();
        // Queue of filtered k-mers.
        let mut q = MonotoneQueue::new();

        // i: position of lmer
        // j: position of kmer
        // j0: position of k0mer

        let o0 = &self.rand_mini.o;

        // 1: init k0-mers.
        let mut k0mers = text.windows(self.k0).enumerate();
        for (j0, k0mer) in k0mers.by_ref().take(self.w0) {
            q0.push(j0, Order::key(o0, k0mer));
        }

        // 2: init k-mers.
        let mut kmers = text.windows(self.k).enumerate().zip(k0mers);
        for ((j, kmer), (j0, k0mer)) in kmers.by_ref().take(self.w - 1) {
            q0.push(j0, Order::key(o0, k0mer));
            let min_pos = q0.pop(j).unwrap().0;
            if min_pos == j || min_pos == j + self.w0 {
                q.push(
                    j,
                    (
                        min_pos == j,
                        Order::key(o0, &text[min_pos..min_pos + self.k0]),
                        self.o.key(kmer),
                    ),
                );
            }
        }

        // 3: Iterate l-mers.
        kmers
            .enumerate()
            .map(move |(i, ((j, kmer), (j0, k0mer)))| {
                q0.push(j0, Order::key(o0, k0mer));
                let min_pos = q0.pop(j).unwrap().0;
                if min_pos == j || min_pos == j + self.w0 {
                    q.push(
                        j,
                        (
                            min_pos == j,
                            Order::key(o0, &text[min_pos..min_pos + self.k0]),
                            self.o.key(kmer),
                        ),
                    );
                }

                q.pop(i).unwrap().0
            })
            .collect()
    }
}

/// Fast %d operations.
/// Taken from https://github.com/lemire/fastmod/blob/master/include/fastmod.h
#[derive(Copy, Clone, Debug)]
pub struct FM32 {
    modulus: u64,
    inv: u64,
}
impl FM32 {
    fn new(modulus: usize) -> Self {
        assert!(modulus <= u32::MAX as usize);
        Self {
            modulus: modulus as u64,
            inv: (u64::MAX / modulus as u64).wrapping_add(1),
        }
    }
    fn reduce(self, h: usize) -> usize {
        // 64-bit multiplication always gives a 128-bit result.
        // All this casting is just needed to extract the high word; it's not slower in practice.
        let lowbits = self.inv.wrapping_mul(h as u64);
        ((lowbits as u128 * self.modulus as u128) >> 64) as usize
    }
}

pub struct ModSampling<O: Order> {
    w: usize,
    l: usize,
    t: usize,
    /// The number of t-mers in an l-mer.
    w_t: usize,
    /// The order on t-mers.
    o: O,
    /// Fast modulo w operations.
    fastmod_w: FM32,
}

impl<O: Order> ModSampling<O> {
    pub fn new(k: usize, w: usize, t: usize, o: O) -> Self {
        assert!(k > 0);
        assert!(t > 0);
        assert!(w > 0);
        let l = k + w - 1;
        let w_t = l - t + 1;
        Self {
            w,
            l,
            t,
            w_t,
            o,
            fastmod_w: FM32::new(w),
        }
    }
}

impl<O: Order> ModSampling<O> {
    pub fn lr_minimizer(k: usize, w: usize, o: O) -> Self {
        assert!(k > w);
        let t = k - w;
        Self::new(k, w, t, o)
    }
    pub fn mod_minimizer(k: usize, w: usize, r: usize, o: O) -> Self {
        let r = r.min(k);
        let t = (k - r) % w + r;
        Self::new(k, w, t, o)
    }
}

impl<O: Order> SamplingScheme for ModSampling<O> {
    fn l(&self) -> usize {
        self.l
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let x = lmer
            .windows(self.t)
            .enumerate()
            .min_by_key(|&(_i, w)| self.o.key(w))
            .unwrap()
            .0;
        x % self.w
    }

    /// NOTE: This is not always a forward scheme.
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        let mut q = MonotoneQueue::new();
        let mut tmers = text.windows(self.t).enumerate();

        for (j, tmer) in tmers.by_ref().take(self.w_t - 1) {
            q.push(j, self.o.key(tmer));
        }
        // i: position of lmer
        // j: position of tmer
        tmers
            .enumerate()
            .map(move |(i, (j, tmer))| {
                q.push(j, self.o.key(tmer));
                i + self.fastmod_w.reduce(q.pop(i).unwrap().0 - i)
            })
            .collect()
    }
}

/// Asymptotic (in k) optimal minimizers:
/// - Assume k=x*w; sum all i mod w positions.
/// - Take minimizer if first coordinate is max or at most sigma-1 away from max.
/// - If multiple, take random.
///
/// Assumes alphabet size sigma=4.
/// Also requires the input to be encoded as values 0..sigma, not as arbitrary ascii codepoints.
///
/// On the overshoot:
/// Instead of the original sigma, we can use sigma-1, since the max difference
/// between two characters is sigma-1.  In fact, I conjecture that
/// floor((sigma+1)/2) overshoot is sufficient. I have a sketch of a proof for
/// this.
pub struct RotMinimizer<O: Order> {
    k: usize,
    w: usize,
    o: O,
}
impl<O: Order> RotMinimizer<O> {
    pub fn new(k: usize, w: usize, o: O) -> Self {
        assert!(k % w == 0);
        Self { k, w, o }
    }
}
impl<O: Order> SamplingScheme for RotMinimizer<O> {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let sigma = 4;
        assert!(lmer.iter().all(|c| *c <= sigma));
        lmer.windows(self.k)
            .enumerate()
            .filter(|(_, kmer)| {
                let psi0 = kmer.iter().step_by(self.w).map(|&x| x as u64).sum::<u64>();
                for j in 1..self.w {
                    let psij = kmer
                        .iter()
                        .skip(j)
                        .step_by(self.w)
                        .map(|&x| x as u64)
                        .sum::<u64>();
                    if !(psij <= psi0 + sigma as u64) {
                        return false;
                    }
                }
                true
            })
            .min_by_key(|&(_i, kmer)| self.o.key(kmer))
            .unwrap()
            .0
    }
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
pub struct AltRotMinimizer<O: Order> {
    k: usize,
    w: usize,
    o: O,
}
impl<O: Order> AltRotMinimizer<O> {
    pub fn new(k: usize, w: usize, o: O) -> Self {
        Self { k, w, o }
    }
}
impl<O: Order> SamplingScheme for AltRotMinimizer<O> {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|&(_i, kmer)| {
                let psi = kmer.iter().step_by(self.w).map(|&x| x as u64).sum::<u64>();
                (Reverse(psi), self.o.key(kmer))
            })
            .unwrap()
            .0
    }
}

pub struct Decycling<O: Order> {
    k: usize,
    w: usize,
    roots: Vec<Complex64>,
    o: O,
    double_decycling: bool,
}

impl<O: Order> Decycling<O> {
    pub fn new(k: usize, w: usize, o: O, double_decycling: bool) -> Self {
        let roots = (0..k)
            .map(|j| Complex64::exp(2. * PI * Complex64::i() / k as f64 * j as f64))
            .collect();
        Self {
            k,
            w,
            roots,
            o,
            double_decycling,
        }
    }
}

impl<O: Order> SamplingScheme for Decycling<O> {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    // TODO: Implement streaming version.
    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|(_i, kmer)| {
                let mut x = Complex64::zero();
                for j in 0..self.k {
                    x += kmer[j] as f64 * self.roots[j];
                }
                let arg = x.arg();
                let key = self.o.key(kmer);
                if PI - 2. * PI / self.k as f64 <= arg && arg < PI {
                    (0, key)
                } else if self.double_decycling && -2. * PI / self.k as f64 <= arg && arg < 0. {
                    (1, key)
                } else {
                    (2, key)
                }
            })
            .unwrap()
            .0
    }
}

/// 1. Find open-syncmer kmers, where the smallest tmer is in the middle.
/// 2. Take the minimizer, preferring open-syncmers.
pub struct OpenSyncmer {
    k: usize,
    w: usize,
    /// Length of tmers to consider inside each kmer.
    t: usize,
    o: RandomOrder,
    rand_mini: Minimizer<RandomOrder>,
    tiebreak: bool,
    closed: bool,
}

impl OpenSyncmer {
    pub fn new(k: usize, w: usize, t: usize, tiebreak: bool, closed: bool) -> Self {
        Self {
            k,
            w,
            t,
            o: RandomOrder,
            rand_mini: Minimizer::new(t, k - t + 1, RandomOrder),
            tiebreak,
            closed,
        }
    }

    fn hash_kmer(&self, kmer: &[u8], j: usize) -> (i32, usize) {
        let m = self.k - self.t;
        let half = m / 2;
        // Prefer kmers whose minimal tmer is in the middle.
        let pref = if j == half {
            0
        } else if self.closed && j == 0 {
            1
        } else if self.closed && j == m {
            1
        } else {
            3
        };
        let tiebreak_hash = if pref == 0 && self.tiebreak {
            let min_tmer = &kmer[j..j + self.t];
            Order::key(&self.rand_mini.o, min_tmer)
        } else {
            Order::key(&self.o, kmer)
        };
        (pref, tiebreak_hash)
    }
}

impl SamplingScheme for OpenSyncmer {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|(_, kmer)| {
                let j = self.rand_mini.sample(kmer);
                self.hash_kmer(kmer, j)
            })
            .unwrap()
            .0
    }

    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        // Queue of t-mers.
        let mut qt = MonotoneQueue::new();
        // Queue of k-mers.
        let mut q = MonotoneQueue::new();

        // i: position of lmer
        // j: position of kmer
        // jt: position of tmer

        let ot = &self.rand_mini.o;

        // 1: init t-mers.
        let mut tmers = text.windows(self.t).enumerate();
        for (jt, tmer) in tmers.by_ref().take(self.k - self.t) {
            qt.push(jt, Order::key(ot, tmer));
        }

        // 2: init k-mers.
        let mut kmers = text.windows(self.k).enumerate().zip(tmers);
        for ((j, kmer), (jt, tmer)) in kmers.by_ref().take(self.w - 1) {
            qt.push(jt, Order::key(ot, tmer));
            q.push(j, self.hash_kmer(kmer, qt.pop(j).unwrap().0 - j));
        }

        // 3: Iterate l-mers.
        kmers
            .enumerate()
            .map(move |(i, ((j, kmer), (jt, tmer)))| {
                qt.push(jt, Order::key(ot, tmer));
                q.push(j, self.hash_kmer(kmer, qt.pop(j).unwrap().0 - j));
                q.pop(i).unwrap().0
            })
            .collect()
    }
}

/// Sample 1/f of the kmers at random, then minimizer to get to 1/w.
pub struct FracMin {
    k: usize,
    w: usize,
    bound: usize,
    o: RandomOrder,
    seed: usize,
}

impl FracMin {
    pub fn new(k: usize, w: usize, f: usize) -> Self {
        let bound = usize::MAX / f;
        Self {
            k,
            w,
            bound,
            o: RandomOrder,
            seed: rand::random(),
        }
    }
}

impl SamplingScheme for FracMin {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|(_, kmer)| {
                let pref = if Order::key(&self.o, kmer) <= self.bound {
                    0
                } else {
                    1
                };
                let kmer_hash = wyhash::wyhash(kmer, self.seed as u64);
                (pref, kmer_hash)
            })
            .unwrap()
            .0
    }
}

/// Prefer first positions t mod w.
/// Then prefer the ends.
pub struct OcModMinimizer<O: Order, OT: Order> {
    k: usize,
    w: usize,
    /// Length of tmers to consider inside each kmer.
    t: usize,
    /// Preferred offset
    offset: usize,
    /// If false, only do the open syncmer part.
    use_closed: bool,
    /// Prefer prefix over suffix?
    prefer_prefix: bool,
    /// Hash middle by tmer.
    open_tmer: bool,
    /// Hash ends by tmer.
    closed_tmer: bool,
    /// Hash remaining cases by tmer
    other_tmer: bool,
    o: O,
    ot: Minimizer<OT>,
}

impl<O: Order, OT: Order> OcModMinimizer<O, OT> {
    pub fn new(
        k: usize,
        w: usize,
        t: usize,
        offset: usize,
        use_closed: bool,
        prefer_prefix: bool,
        open_tmer: bool,
        closed_tmer: bool,
        other_tmer: bool,
        o: O,
        ot: OT,
    ) -> Self {
        Self {
            k,
            w,
            t,
            offset,
            use_closed,
            prefer_prefix,
            open_tmer,
            closed_tmer,
            other_tmer,
            o,
            ot: Minimizer::new(t, k - t + 1, ot),
        }
    }

    fn hash_kmer(&self, kmer: &[u8], j: usize) -> (i32, usize) {
        assert!(j <= self.k - self.t);
        // Prefer kmers whose minimal tmer is in the middle.
        let pref = if j % self.w == self.offset % self.w {
            0
        } else if self.use_closed && j == 0 {
            1
        } else if self.use_closed && j == self.k - self.t {
            if self.prefer_prefix {
                2
            } else {
                1
            }
        } else {
            3
        };
        let hash_tmer = if pref == 0 {
            self.open_tmer
        } else if pref == 1 || pref == 2 {
            self.closed_tmer
        } else {
            self.other_tmer
        };
        let tiebreak_hash = if hash_tmer {
            let min_tmer = &kmer[j..j + self.t];
            Order::key(&self.ot.o, min_tmer)
        } else {
            Order::key(&self.o, kmer)
        };
        (pref, tiebreak_hash)
    }
}

impl<O: Order, OT: Order> SamplingScheme for OcModMinimizer<O, OT> {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|(_, kmer)| {
                let j = self.ot.sample(kmer);
                self.hash_kmer(kmer, j)
            })
            .unwrap()
            .0
    }

    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        // Queue of t-mers.
        let mut qt = MonotoneQueue::new();
        // Queue of k-mers.
        let mut q = MonotoneQueue::new();

        // i: position of lmer
        // j: position of kmer
        // jt: position of tmer

        let ot = &self.ot.o;

        // 1: init t-mers.
        let mut tmers = text.windows(self.t).enumerate();
        for (jt, tmer) in tmers.by_ref().take(self.k - self.t) {
            qt.push(jt, Order::key(ot, tmer));
        }

        // 2: init k-mers.
        let mut kmers = text.windows(self.k).enumerate().zip(tmers);
        for ((j, kmer), (jt, tmer)) in kmers.by_ref().take(self.w - 1) {
            qt.push(jt, Order::key(ot, tmer));
            q.push(j, self.hash_kmer(kmer, qt.pop(j).unwrap().0 - j));
        }

        // 3: Iterate l-mers.
        kmers
            .enumerate()
            .map(move |(i, ((j, kmer), (jt, tmer)))| {
                qt.push(jt, Order::key(ot, tmer));
                q.push(j, self.hash_kmer(kmer, qt.pop(j).unwrap().0 - j));
                q.pop(i).unwrap().0
            })
            .collect()
    }
}

pub struct ClosedSyncmer {
    k: usize,
    w: usize,
    /// Length of tmers to consider inside each kmer.
    t: usize,
    h: usize,
    loose: bool,
    open: bool,
    o: RandomOrder,
    rand_mini: Minimizer<RandomOrder>,
    threshold: usize,
    threshold2: usize,
    count: Cell<usize>,
    count2: Cell<usize>,
    count3: Cell<usize>,
}

impl ClosedSyncmer {
    pub fn new(k: usize, w: usize, t: usize, loose: bool, open: bool, h: usize) -> Self {
        // let k2 = k - t + 1;
        let threshold = t * h;
        let threshold2 = (2 * t - 256) * h;
        eprintln!("Thresholds: {} {}", threshold, threshold2);
        Self {
            k,
            w,
            t,
            h,
            loose,
            open,
            o: RandomOrder,
            // rand_mini: Minimizer::new(t, k - t + 1, RandomOrder),
            rand_mini: Minimizer::new(1, k, RandomOrder),
            threshold,
            threshold2,
            count: Cell::new(0),
            count2: Cell::new(0),
            count3: Cell::new(0),
        }
    }

    fn hash_kmer(&self, kmer: &[u8]) -> (i32, isize, usize) {
        // Hash all tmers in the kmer.
        // If the first tmer is the only one > threshold, it's a good kmer.
        let mut pref = 0;
        // if Order::key(&self.rand_mini.o, &kmer[0..self.t]) >= self.threshold {
        //     for tmer in kmer[1..].windows(self.t) {
        //         if Order::key(&self.rand_mini.o, tmer) < self.threshold {
        //             pref = 1;
        //             break;
        //         }
        //     }
        // } else {
        //     pref = 1;
        // }
        pref = 100;
        let e = self.h.min(self.k);

        let s0: usize = kmer[0..e].iter().map(|&c| c as usize).sum();
        if s0 >= self.threshold {
            pref -= 1;
        }

        if pref < 100 {
            for hmer in kmer[1..].windows(self.h) {
                if hmer.iter().map(|&c| c as usize).sum::<usize>() < self.threshold {
                    pref -= 1;
                } else {
                    break;
                }
            }
        }

        if false {
            // Right-closed kmers
            // if pref > 100 - self.k as i32 {
            //     pref = 200;
            //     if (kmer[self.k - 1] as usize) >= self.threshold {
            //         pref -= 1;
            //     }

            //     if pref < 200 {
            //         for hmer in kmer[..self.k - 1].windows(self.h).rev() {
            //             if hmer.iter().map(|&c| c as usize).sum::<usize>() < self.threshold {
            //                 pref -= 1;
            //             } else {
            //                 break;
            //             }
            //         }
            //     }

            //     if !self.loose && pref > 200 - self.k as i32 {
            //         pref = 200;
            //     }
            // }
        } else {
            if !self.loose && pref > 100 - self.k as i32 {
                pref = 100;
            }
        }

        // if self.open && pref == 100 {
        //     pref = 200;
        //     let s0: usize = kmer[0..e].iter().map(|&c| c as usize).sum();
        //     if s0 >= self.threshold2 {
        //         pref -= 1;
        //     }

        //     if pref < 200 {
        //         for hmer in kmer[1..].windows(self.h) {
        //             if hmer.iter().map(|&c| c as usize).sum::<usize>() < self.threshold2 {
        //                 pref -= 1;
        //             } else {
        //                 break;
        //             }
        //         }
        //     }
        // }

        // if pref == 1 {
        // right-max kmer?
        // if (kmer[self.k - 1] as usize) < self.threshold {
        //     pref = 2;
        // }
        // for &tmer in &kmer[..self.k - 1] {
        //     if (tmer as usize) >= self.threshold {
        //         pref = 2;
        //     }
        // }
        // }
        // let v = kmer[0] as isize - *kmer[1..].iter().max().unwrap() as isize;
        // } else {
        //     pref = 1;
        // }
        //
        // Classic max-in-middle open variant.
        //
        // if self.open && pref as usize >= 100 {
        //     let m = self.k / 2;
        //     for (i, &tmer) in kmer.iter().enumerate() {
        //         if ((tmer as usize) >= self.threshold) ^ (i == m) {
        //             pref = 200;
        //         }
        //     }
        // }
        (
            pref,
            0,
            // -v,
            // (255 - kmer[0]) as isize,
            // kmer[0] as isize - m as isize,
            Order::key(&self.o, kmer),
        )
    }
}

impl SamplingScheme for ClosedSyncmer {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let mut good = false;

        let h = self.hash_kmer(&lmer[0..self.k]);
        if h.0 < 100 {
            self.count3.set(self.count3.get() + 1);
        }

        // Do not allow kmers that are followed by only 0s, since they are likely to overlap a 10000.
        let mut skip = 0;
        if self.open {}

        let ((pref, _, _), mut p, _) = lmer
            .windows(self.k)
            .enumerate()
            .map(|(i, kmer)| {
                let h = self.hash_kmer(kmer);
                if h.0 < 100 {
                    self.count.set(self.count.get() + 1);
                    good = true;
                }
                (h, i, kmer)
            })
            .min()
            .unwrap();

        // if self.open && pref == 100 && self.k <= self.w {
        //     // Check if maybe there is a large char in the last k-1 positions.
        //     // Then take the kmer ending just before it.
        //     for (i, &c) in lmer[self.w..].iter().enumerate() {
        //         if c >= self.threshold as u8 {
        //             p = i + self.w - self.k;
        //         }
        //     }
        // }

        if good {
            self.count2.set(self.count2.get() + 1);
        }
        p
    }
}

impl Drop for ClosedSyncmer {
    fn drop(&mut self) {
        eprintln!("k {} => count : {}", self.k, self.count.get());
        eprintln!("k {} => count2: {}", self.k, self.count2.get());
        eprintln!("k {} => count3: {}", self.k, self.count3.get());
    }
}

// k = 7
// ---*---
//      ---*---
// 52148367
// --*--
//
// 7163524
//   ^^^   w-k
//
// 1.24.59.b.a8.76.3.
// --*--

// w=4
// k=5= w+1
// l=8=2w
// c=9=2w+1
//   016238475
//   *   *   *
//    ... ...
//     0
//     *----
//   density of (w, w+1) <= 3/(2w+1) ~ 3/(2w) = 1.5/w
//
//   1.67/w
//
//   016238475
//   *----     charge
//       ----* charge
//   ----*     charge
//       *---- charge
//
// 0 *----
// 1  *----
// -   -*---
// 2    *----
// 3     *----
