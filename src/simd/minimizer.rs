use crate::{simd::linearize, Captures};

use super::{
    collect::{collect, collect_and_dedup},
    dedup,
    nthash::{nthash32_par_it, nthash32_scalar_it},
};
use core::array::from_fn;
use itertools::Itertools;
use packed_seq::{Seq, S};

/// A custom RingBuf implementation that has a fixed size `w` and wraps around.
struct RingBuf<V> {
    w: usize,
    idx: usize,
    data: Vec<V>,
}

impl<V: Clone> RingBuf<V> {
    #[inline(always)]
    fn new(w: usize, v: V) -> Self {
        assert!(w > 0);
        let data = vec![v; w];
        RingBuf { w, idx: 0, data }
    }

    /// Returns the next index to be written.
    #[inline(always)]
    fn idx(&self) -> usize {
        self.idx
    }

    #[inline(always)]
    fn push(&mut self, v: V) {
        self.data[self.idx] = v;
        self.idx += 1;
        if self.idx == self.w {
            self.idx = 0;
        }
    }
}

/// A RingBuf can be used as a slice.
impl<V> std::ops::Deref for RingBuf<V> {
    type Target = [V];

    #[inline(always)]
    fn deref(&self) -> &[V] {
        &self.data
    }
}

/// A RingBuf can be used as a mutable slice.
impl<V> std::ops::DerefMut for RingBuf<V> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut [V] {
        &mut self.data
    }
}

#[inline(always)]
fn sliding_min_scalar_it(
    it: impl ExactSizeIterator<Item = u32>,
    w: usize,
) -> impl ExactSizeIterator<Item = u32> {
    assert!(w > 0);
    assert!(w < (1 << 15), "This method is not tested for large w.");
    assert!(it.size_hint().0 < (1 << 32));
    let mut prefix_min = u32::MAX;
    let mut ring_buf = RingBuf::new(w, prefix_min);
    // We only compare the upper 16 bits of each hash.
    // Ties are broken automatically in favour of lower pos.
    let val_mask = 0xffff_0000;
    let pos_mask = 0x0000_ffff;
    let mut pos = 0;
    let max_pos = (1 << 16) - 1;
    let mut pos_offset = 0;

    let mut it = it.map(
        #[inline(always)]
        move |val| {
            // Make sure the position does not interfere with the hash value.
            if pos == max_pos {
                let delta = ((1 << 16) - 2 - w) as u32;
                pos -= delta;
                prefix_min -= delta;
                pos_offset += delta;
                for x in &mut *ring_buf {
                    *x -= delta;
                }
            }
            let elem = (val & val_mask) | pos;
            pos += 1;
            ring_buf.push(elem);
            prefix_min = prefix_min.min(elem);
            // After a chunk has been filled, compute suffix minima.
            if ring_buf.idx() == 0 {
                let mut suffix_min = ring_buf[w - 1];
                for i in (0..w - 1).rev() {
                    suffix_min = suffix_min.min(ring_buf[i]);
                    ring_buf[i] = suffix_min;
                }
                prefix_min = elem; // slightly faster than assigning S::splat(u32::MAX)
            }
            let suffix_min = unsafe { *ring_buf.get_unchecked(ring_buf.idx()) };
            (prefix_min.min(suffix_min) & pos_mask) + pos_offset
        },
    );
    // This optimizes better than it.skip(w-1).
    it.by_ref().take(w - 1).for_each(drop);
    it
}

#[inline(always)]
fn sliding_min_par_it(
    it: impl ExactSizeIterator<Item = S>,
    w: usize,
) -> impl ExactSizeIterator<Item = S> {
    assert!(w > 0);
    assert!(w < (1 << 15), "This method is not tested for large w.");
    assert!(it.size_hint().0 * 8 < (1 << 32));
    let mut prefix_min = S::splat(u32::MAX);
    let mut ring_buf = RingBuf::new(w, prefix_min);
    // We only compare the upper 16 bits of each hash.
    // Ties are broken automatically in favour of lower pos.
    let val_mask = S::splat(0xffff_0000);
    let pos_mask = S::splat(0x0000_ffff);
    let max_pos = S::splat((1 << 16) - 1);
    let mut pos = S::splat(0);
    let mut pos_offset: S =
        from_fn(|l| (l * (it.size_hint().0.saturating_sub(w - 1))) as u32).into();

    let mut it = it.map(
        #[inline(always)]
        move |val| {
            // Make sure the position does not interfere with the hash value.
            if pos == max_pos {
                let delta = S::splat((1 << 16) - 2 - w as u32);
                pos -= delta;
                prefix_min -= delta;
                pos_offset += delta;
                for x in &mut *ring_buf {
                    *x -= delta;
                }
            }
            let elem = (val & val_mask) | pos;
            pos += S::splat(1);
            ring_buf.push(elem);
            prefix_min = prefix_min.min(elem);
            // After a chunk has been filled, compute suffix minima.
            if ring_buf.idx() == 0 {
                let mut suffix_min = ring_buf[w - 1];
                for i in (0..w - 1).rev() {
                    suffix_min = suffix_min.min(ring_buf[i]);
                    ring_buf[i] = suffix_min;
                }
                prefix_min = elem; // slightly faster than assigning S::splat(u32::MAX)
            }
            let suffix_min = unsafe { *ring_buf.get_unchecked(ring_buf.idx()) };
            (prefix_min.min(suffix_min) & pos_mask) + pos_offset
        },
    );
    // This optimizes better than it.skip(w-1).
    it.by_ref().take(w - 1).for_each(drop);
    it
}

/// Returns the minimizer of a window using a naive linear scan.
/// Uses NT hash with canonical hashes when `RC` is true.
pub fn minimizer_window_naive<'s, const RC: bool>(seq: impl Seq<'s>, k: usize) -> usize {
    nthash32_scalar_it::<RC>(seq, k)
        .map(|x| x & 0xffff_0000)
        .position_min()
        .unwrap()
}

/// Returns an iterator over the absolute positions of the minimizers of a sequence.
/// Returns one value for each window of size `w+k-1` in the input. Use
/// `Itertools::dedup()` to obtain the distinct positions of the minimizers.
///
/// Prefer `minimizer_simd_it` that internally used SIMD, or `minimizer_par_it` if it works for you.
pub fn minimizer_scalar_it<'s, const RC: bool>(
    seq: impl Seq<'s>,
    k: usize,
    w: usize,
) -> impl ExactSizeIterator<Item = u32> + Captures<&'s ()> {
    let it = nthash32_scalar_it::<RC>(seq, k);
    sliding_min_scalar_it(it, w)
}

/// Returns an iterator over the absolute positions of the minimizers of a sequence.
/// Returns one value for each window of size `w+k-1` in the input. Use
/// `Itertools::dedup()` to obtain the distinct positions of the minimizers.
///
/// This splits the windows of the sequence into chunks of 2^13.
/// Minimizers for each chunk are eagerly computed using 8 parallel streams using SIMD using `minimizers_par_it`.
/// Then returns a linear iterator over the buffer.
/// Once the buffer runs out, the next chunk is computed.
///
/// NOTE: This method is ~4x slower than the minimizer computation itself, and
///       only ~2x faster than the scalar version. Mostly because shuffling memory is slow.
/// TODO: Fix this.
pub fn minimizer_simd_it<'s, const RC: bool>(
    seq: impl Seq<'s>,
    k: usize,
    w: usize,
) -> impl ExactSizeIterator<Item = u32> + Captures<&'s ()> {
    linearize::linearize_with_offset(seq, k + w - 1, move |seq| minimizer_par_it::<RC>(seq, k, w))
}

/// Split the windows of the sequence into 8 chunks of equal length ~len/8.
/// Then return the positions of the minimizers of each of them in parallel using SIMD,
/// and return the remaining few using the second iterator.
// TODO: Take a hash function as argument.
pub fn minimizer_par_it<'s, const RC: bool>(
    seq: impl Seq<'s>,
    k: usize,
    w: usize,
) -> (
    impl ExactSizeIterator<Item = S> + Captures<&'s ()>,
    impl ExactSizeIterator<Item = u32> + Captures<&'s ()>,
) {
    let (par_head, tail) = nthash32_par_it::<RC>(seq, k, w);
    let par_head = sliding_min_par_it(par_head, w);
    let offset = 8 * par_head.size_hint().0 as u32;
    let tail = sliding_min_scalar_it(tail, w).map(move |pos| offset + pos);
    (par_head, tail)
}

pub fn minimizers_collect<'s, const RC: bool>(seq: impl Seq<'s>, k: usize, w: usize) -> Vec<u32> {
    let (par_head, tail) = minimizer_par_it::<RC>(seq, k, w);
    collect(par_head, tail)
}

pub fn minimizers_dedup<'s, const RC: bool>(seq: impl Seq<'s>, k: usize, w: usize) -> Vec<u32> {
    let (par_head, tail) = minimizer_par_it::<RC>(seq, k, w);
    let mut positions = collect(par_head, tail);
    dedup(&mut positions);
    positions
}

pub fn minimizers_collect_and_dedup<'s, const RC: bool>(
    seq: impl Seq<'s>,
    k: usize,
    w: usize,
) -> Vec<u32> {
    let (par_head, tail) = minimizer_par_it::<RC>(seq, k, w);
    collect_and_dedup(par_head, tail)
}

#[cfg(test)]
mod test {
    use packed_seq::{AsciiSeq, AsciiSeqVec, PackedSeqVec, SeqVec};

    use super::*;
    use std::{cell::LazyCell, iter::once};

    const ASCII_SEQ: LazyCell<AsciiSeqVec> = LazyCell::new(|| AsciiSeqVec::random(1024 * 1024));
    const PACKED_SEQ: LazyCell<PackedSeqVec> = LazyCell::new(|| PackedSeqVec::random(1024 * 1024));

    #[test]
    fn scalar_byte() {
        let seq = &*ASCII_SEQ;
        for k in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
            for w in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
                for len in (0..100).chain(once(1024 * 32)) {
                    let seq = seq.slice(0..len);
                    let single = seq.0[0..len]
                        .windows(w + k - 1)
                        .enumerate()
                        .map(|(pos, seq)| {
                            (pos + minimizer_window_naive::<false>(AsciiSeq::new(seq, k), k)) as u32
                        })
                        .collect::<Vec<_>>();
                    let scalar = minimizer_scalar_it::<false>(seq, k, w).collect::<Vec<_>>();
                    assert_eq!(single, scalar, "k={k}, w={w}, len={len}");
                }
            }
        }
    }

    #[test]
    fn parallel_iter_byte() {
        let seq = &*ASCII_SEQ;
        for k in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
            for w in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
                for len in (0..100).chain(once(1024 * 128)) {
                    let seq = seq.slice(0..len);
                    let scalar = minimizer_scalar_it::<false>(seq, k, w).collect::<Vec<_>>();
                    let (par_head, tail) = minimizer_par_it::<false>(seq, k, w);
                    let par_head = par_head.collect::<Vec<_>>();
                    let parallel_iter = (0..8)
                        .flat_map(|l| par_head.iter().map(move |x| x.as_array_ref()[l]))
                        .chain(tail)
                        .collect::<Vec<_>>();

                    assert_eq!(scalar, parallel_iter, "k={k}, w={w}, len={len}");
                }
            }
        }
    }

    #[test]
    fn parallel_iter_packed() {
        let seq = &*PACKED_SEQ;
        for k in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
            for w in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
                for len in (0..100).chain(once(1024 * 128)) {
                    let seq = seq.slice(0..len);
                    let scalar = minimizer_scalar_it::<false>(seq, k, w).collect::<Vec<_>>();
                    let (par_head, tail) = minimizer_par_it::<false>(seq, k, w);
                    let par_head = par_head.collect::<Vec<_>>();
                    let parallel_iter = (0..8)
                        .flat_map(|l| par_head.iter().map(move |x| x.as_array_ref()[l]))
                        .chain(tail)
                        .collect::<Vec<_>>();
                    assert_eq!(scalar, parallel_iter, "k={k}, w={w}, len={len}");
                }
            }
        }
    }

    #[test]
    fn linearized() {
        let seq = &*PACKED_SEQ;
        for k in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
            for w in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
                for len in (0..100).chain(once(1024 * 128 + 765)) {
                    let seq = seq.slice(0..len);
                    let scalar = minimizer_scalar_it::<false>(seq, k, w).collect::<Vec<_>>();
                    let simd = minimizer_simd_it::<false>(seq, k, w).collect::<Vec<_>>();
                    assert_eq!(
                        scalar,
                        simd,
                        "k={}, len={} len left {} len right {}",
                        k,
                        len,
                        scalar.len(),
                        simd.len()
                    );
                }
            }
        }
    }
}
