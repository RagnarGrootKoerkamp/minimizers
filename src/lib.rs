// TODO: Implement fast version of AntiLex sort, skipping NtHash.
// TODO: Make the 'linearized' `minimizers_simd_it` and `nthash32_simd_it` iterators faster.
//       Currently they are 4x slower than the `*_par_it` versions.

pub mod de_bruijn_seq;
pub mod monotone_queue;
pub mod order;
#[cfg(feature = "python-bindings")]
mod py;
pub mod schemes;
pub mod simd;

use crate::monotone_queue::MonotoneQueue;
use itertools::Itertools;
use order::*;
use rand_chacha::{
    rand_core::{RngCore, SeedableRng},
    ChaChaRng,
};
use std::fmt::Debug;

/// An iterator over *all* minimizer positions. Not yet deduplicated.
///
/// NOTE: For non-forward schemes, positions may be returned twice.
pub trait MinimizerIt: Iterator<Item = usize> {}
impl<I: Iterator<Item = usize>> MinimizerIt for I {}

#[typetag::serialize(tag = "minimizer_type")]
pub trait Params: Debug + Sync {
    fn build(&self, w: usize, k: usize, sigma: usize) -> Box<dyn SamplingScheme>;
}

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

pub trait Captures<U> {}
impl<T: ?Sized, U> Captures<U> for T {}
