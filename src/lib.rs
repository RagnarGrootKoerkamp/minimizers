#![feature(type_alias_impl_trait, iter_repeat_n, trait_alias)]

pub mod bruteforce;
pub mod de_bruijn_seq;
pub mod ilp_scheme;
pub mod local_scheme;
pub mod monotone_queue;
pub mod order;

use crate::monotone_queue::MonotoneQueue;
use itertools::Itertools;
use num::{complex::Complex64, Zero};
use order::*;
use std::{cmp::Reverse, f64::consts::PI};

/// An iterator over deduplicated minimizer positions.
///
/// NOTE: For non-forward schemes, positions may be returned twice.
pub trait MinimizerIt = Iterator<Item = usize>;

pub trait SamplingScheme {
    fn l(&self) -> usize;
    /// Sample a single lmer.
    fn sample(&self, lmer: &[u8]) -> usize;
    /// Sample all lmers in a text, and deduplicate.
    /// This default implementation simply calls `sample` on each lmer.
    // TODO: Take an iterator over u8 instead?
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        self.stream_naive(text)
    }
    /// Sample all lmers in a cyclic text of length `len`.
    /// Text must have length at least `len + l-1`, and the additional
    /// characters must equal a prefix of the text.
    ///
    /// Returns the number of selected minimizers.
    fn sample_cyclic_text(&self, text: &[u8], len: usize) -> usize {
        let text = &text[..len + self.l() - 1];
        let mut poss: Vec<_> = self.stream(text).map(|p| p % len).collect();
        poss.sort();
        poss.dedup();
        poss.len()
    }
    /// Sample all lmers in a text, and deduplicate.
    /// This default implementation simply calls `sample` on each lmer.
    fn stream_naive(&self, text: &[u8]) -> impl MinimizerIt {
        text.windows(self.l())
            .enumerate()
            .map(|(i, lmer)| i + self.sample(lmer))
            .dedup()
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
        Self {
            k,
            w,
            l: k + w - 1,
            o,
        }
    }
    pub fn ord(&self) -> &O {
        &self.o
    }
}

impl<O: DirectedOrder> SamplingScheme for Minimizer<O> {
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
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
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
            .dedup()
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
    l: usize,
    r: usize,
}

impl BdAnchor {
    pub fn new(l: usize, r: usize) -> Self {
        Self { l, r }
    }
}

impl SamplingScheme for BdAnchor {
    fn l(&self) -> usize {
        self.l
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        debug_assert_eq!(lmer.len(), self.l());
        let mut best = 0;
        for i in 1..lmer.len().saturating_sub(self.r) {
            for j in 0..lmer.len() {
                if lmer[(i + j) % lmer.len()] != lmer[(best + j) % lmer.len()] {
                    if lmer[(i + j) % lmer.len()] < lmer[(best + j) % lmer.len()] {
                        best = i;
                    }
                    break;
                }
            }
        }
        best
    }
}

pub struct Miniception<O: Order> {
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

impl<O: Order> Miniception<O> {
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

impl<O: Order> SamplingScheme for Miniception<O> {
    fn l(&self) -> usize {
        self.l
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .filter(|(_, kmer)| {
                let j = self.rand_mini.sample(kmer);
                j == 0 || j == self.k - self.k0
            })
            .min_by_key(|&(_i, kmer)| self.o.key(kmer))
            .unwrap()
            .0
    }

    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
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
            .dedup()
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

    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
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
            .dedup()
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

impl ModSampling<RandomOrder> {
    pub fn lr_minimizer(k: usize, w: usize) -> Self {
        assert!(k > w);
        let t = k - w;
        Self::new(k, w, t, RandomOrder)
    }
    pub fn mod_minimizer(k: usize, w: usize) -> Self {
        let r = 5;
        assert!(k > r);
        let t = (k - r) % w + r;
        Self::new(k, w, t, RandomOrder)
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
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
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
            .dedup()
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
