#![feature(exclusive_range_pattern, type_alias_impl_trait)]

pub mod monotone_queue;

use crate::monotone_queue::MonotoneQueue;
use itertools::Itertools;
use std::{f32::consts::PI, iter::zip};

pub trait Minimizers: Iterator<Item = usize> {}
impl<T: Iterator<Item = usize>> Minimizers for T {}

fn h(kmer: &[u8]) -> u64 {
    fxhash::hash64(kmer)
}

// TODO: Make `random_minimizer` take an `Order`.

pub fn random_minimizer(s: &[u8], k: usize) -> usize {
    assert!(k > 0);
    assert!(k <= s.len());
    s.windows(k)
        .enumerate()
        .min_by_key(|&(_i, kmer)| h(kmer))
        .unwrap()
        .0
}

pub fn text_minimizers<'a>(text: &'a [u8], w: usize, k: usize) -> impl Minimizers + 'a {
    let mut q = monotone_queue::MonotoneQueue::new();
    let mut kmers = text.windows(k).enumerate();
    for (j, kmer) in kmers.by_ref().take(w - 1) {
        q.push(j, h(kmer));
    }
    // i: position of lmer
    // j: position of kmer
    kmers.enumerate().map(move |(i, (j, kmer))| {
        q.push(j, h(kmer));
        q.pop(i).unwrap().0
    })
}

pub fn bd_anchor(s: &[u8], r: usize) -> usize {
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
pub fn miniception(s: &[u8], k: usize, k0: usize) -> usize {
    let w = s.len() - k + 1;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    let w0 = k - k0;
    let _l0 = w0 + k0 - 1;

    s.windows(k)
        .enumerate()
        .filter(|(_, kmer)| {
            let i = random_minimizer(kmer, k0);
            assert!(i <= k - k0);
            i == 0 || i == k - k0
        })
        .min_by_key(|&(_i, w)| h(w))
        .unwrap()
        .0
}

pub fn text_miniception<'a>(text: &'a [u8], w: usize, k: usize, k0: usize) -> impl Minimizers + 'a {
    // The number of k0-mers in a kmer.
    let w0 = k - k0;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    // Queue of all k0-mers.
    let mut q0 = MonotoneQueue::new();
    // Queue of filtered k-mers.
    let mut q = MonotoneQueue::new();

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
    kmers.enumerate().map(move |(i, ((j, kmer), (j0, k0mer)))| {
        q0.push(j0, h(k0mer));
        let min_pos = q0.pop(j).unwrap().0;
        if min_pos == j || min_pos == j + w0 {
            q.push(j, h(kmer));
        }

        q.pop(i).unwrap().0
    })
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
pub fn miniception_new(s: &[u8], k: usize, k0: usize) -> usize {
    let w = s.len() - k + 1;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    let w0 = k - k0;
    let _l0 = w0 + k0 - 1;

    s.windows(k)
        .enumerate()
        .filter_map(|(i, kmer)| {
            let j = random_minimizer(kmer, k0);
            assert!(j <= k - k0);
            if j == 0 || j == k - k0 {
                Some((i, kmer, j == 0, h(&kmer[j..j + k0])))
            } else {
                None
            }
        })
        .min_by_key(|&(_i, w, b, hk0)| (b, hk0, h(w)))
        .unwrap()
        .0
}

// TODO: dedup
pub fn text_miniception_new<'a>(
    text: &'a [u8],
    w: usize,
    k: usize,
    k0: usize,
) -> impl Minimizers + 'a {
    // The number of k0-mers in a kmer.
    let w0 = k - k0;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    // Queue of all k0-mers.
    let mut q0 = MonotoneQueue::new();
    // Queue of filtered k-mers.
    let mut q = MonotoneQueue::new();

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
    kmers.enumerate().map(move |(i, ((j, kmer), (j0, k0mer)))| {
        q0.push(j0, h(k0mer));
        let min_pos = q0.pop(j).unwrap().0;
        if min_pos == j || min_pos == j + w0 {
            q.push(j, (min_pos == j, h(&text[min_pos..min_pos + k0]), h(kmer)));
        }

        q.pop(i).unwrap().0
    })
}

pub fn robust_biminimizer(s: &[u8], k: usize, last: &mut usize) -> usize {
    let mut vals = s
        .windows(k)
        .enumerate()
        .map(|(i, w)| (h(w), i))
        .collect_vec();
    vals.sort();
    let i1 = vals[0].1;
    if vals.len() == 1 {
        *last = i1;
        return i1;
    }
    let i2 = vals[1].1;
    if *last == i1 + 1 || *last == i2 + 1 {
        *last -= 1;
    } else {
        *last = i1.max(i2);
    }
    *last
}

/// Find minimal t-mer at pos idx. Then select idx % w.
/// NOTE: THIS IS NOT A FORWARD SCHEME.
pub fn lr_minimizer(s: &[u8], k: usize, t: usize) -> usize {
    let l = s.len();
    let w = l - k + 1;
    let idx = s
        .windows(t)
        .enumerate()
        .min_by_key(|&(_i, w)| h(w))
        .unwrap()
        .0;
    let i = if idx >= w { idx - w } else { idx };
    assert!(i + k <= l);
    i
}

/// NOTE: THIS IS NOT A FORWARD SCHEME.
pub fn text_lr_minimizers<'a>(
    text: &'a [u8],
    w: usize,
    k: usize,
    t: usize,
) -> impl Minimizers + 'a {
    let mut q = MonotoneQueue::new();
    let l = w + k - 1;
    let wt = l - t + 1;
    let mut tmers = text.windows(t).enumerate();

    for (j, tmer) in tmers.by_ref().take(wt - 1) {
        q.push(j, h(tmer));
    }
    // i: position of lmer
    // j: position of tmer
    tmers.enumerate().map(move |(i, (j, tmer))| {
        q.push(j, h(tmer));
        let idx = q.pop(i).unwrap().0 - i;
        i + if idx >= w { idx - w } else { idx }
    })
}

/// Find minimal t-mer at pos idx. Then select idx % w.
/// NOTE: THIS IS NOT A FORWARD SCHEME.
pub fn mod_minimizer(s: &[u8], k: usize, t: usize) -> usize {
    let l = s.len();
    let w = l - k + 1;
    let idx = s
        .windows(t)
        .enumerate()
        .min_by_key(|&(_i, w)| h(w))
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
            m: (u64::MAX / d as u64).wrapping_add(1),
        }
    }
    fn reduce(self, h: usize) -> usize {
        let lowbits = self.m.wrapping_mul(h as u64);
        ((lowbits as u128 * self.d as u128) >> 64) as usize
    }
}

/// NOTE: THIS IS NOT A FORWARD SCHEME.
pub fn text_mod_minimizers<'a>(
    text: &'a [u8],
    w: usize,
    k: usize,
    t: usize,
) -> impl Minimizers + 'a {
    let mut q = MonotoneQueue::new();
    let l = w + k - 1;
    let wt = l - t + 1;
    let mut tmers = text.windows(t).enumerate();

    let fastmod_w = FM32::new(w);

    for (j, tmer) in tmers.by_ref().take(wt - 1) {
        q.push(j, h(tmer));
    }
    // i: position of lmer
    // j: position of tmer
    tmers.enumerate().map(move |(i, (j, tmer))| {
        q.push(j, h(tmer));
        i + fastmod_w.reduce(q.pop(i).unwrap().0 - i)
    })
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
pub fn rot_minimizer(s: &[u8], k: usize) -> usize {
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

pub fn decycling_minimizer_init(k: usize) -> Vec<f32> {
    // (0..k).map(|i| (i as f32).sin()).collect()
    (0..k)
        .map(|i| (i as f32 / k as f32 * 2. * PI).sin())
        .collect()
}

/// Check if this is in the decyling set.
/// Note that we ignore the edge case where ix=0 and ix2=0.
pub fn is_in_decycling_set(kmer: &[u8], cs: &Vec<f32>) -> bool {
    let ix = zip(kmer, cs).map(|(&x, c)| x as f32 * c).sum::<f32>();
    let ix2 = *kmer.last().unwrap() as f32 * cs[0]
        + zip(kmer, &cs[1..]).map(|(&x, c)| x as f32 * c).sum::<f32>();
    ix > 0. && ix2 <= 0.
}
pub fn is_in_double_decycling_set(kmer: &[u8], cs: &Vec<f32>) -> bool {
    let ix = zip(kmer, cs).map(|(&x, c)| x as f32 * c).sum::<f32>();
    let ix2 = *kmer.last().unwrap() as f32 * cs[0]
        + zip(kmer, &cs[1..]).map(|(&x, c)| x as f32 * c).sum::<f32>();
    ix < 0. && ix2 >= 0.
}

/// Decycling minimizer using Mykkelvelt embedding
pub fn decycling_minimizer(s: &[u8], k: usize, cs: &Vec<f32>) -> usize {
    s.windows(k)
        .enumerate()
        .min_by_key(|&(_i, kmer)| (!is_in_decycling_set(kmer, cs), h(kmer)))
        .unwrap()
        .0
}

/// Decycling minimizer using Mykkelvelt embedding
pub fn double_decycling_minimizer(s: &[u8], k: usize, cs: &Vec<f32>) -> usize {
    s.windows(k)
        .enumerate()
        .min_by_key(|&(_i, kmer)| {
            (
                !is_in_decycling_set(kmer, cs),
                !is_in_double_decycling_set(kmer, cs),
                h(kmer),
            )
        })
        .unwrap()
        .0
}

pub struct SuperKmer {
    /// Start position in text of the first basepair.
    pub start: usize,
    /// Past-the-end index of the start of the last kmer.
    pub kmer_end: usize,
    /// Past-the-end index of the full super-kmer.
    pub super_kmer_end: usize,
    /// The position of the minimizer in the super-kmer.
    pub minimizer_pos: usize,
}

/// Takes an iterator over the minimize positions in a string and returns an
/// iterator over super kmers, i.e. substrings sharing the same minimizer.
///
/// Returns (start, end, &[u8])
pub fn super_kmers<'a>(
    text: &'a [u8],
    k: usize,
    minimizer_positions: impl Minimizers + 'a,
) -> impl Iterator<Item = SuperKmer> + 'a {
    let mut first = None;
    minimizer_positions.enumerate().filter_map(move |(i, pos)| {
        let Some((first_i, first_pos)) = first else {
            // Initial set.
            first = Some((i, pos));
            return None;
        };
        if i + 1 == text.len() {
            // Save last.
            return Some(SuperKmer {
                start: first_i,
                kmer_end: i + 1,
                super_kmer_end: i + 1 + k,
                minimizer_pos: first_pos,
            });
        }
        // Extend?
        if pos == first_pos {
            return None;
        } else {
            // Start new.
            first = Some((i, pos));
            // Save previous.
            return Some(SuperKmer {
                start: first_i,
                kmer_end: i,
                super_kmer_end: i + k,
                minimizer_pos: first_pos,
            });
        }
    })
}
