use std::{
    array::from_fn,
    mem::transmute,
    simd::{ptr::SimdConstPtr, Simd},
};

use num::Integer;

pub(crate) const L: usize = 8;
pub(crate) type S = Simd<u32, L>;

/// Iterate over the 2-bit characters of a sequence type.
pub trait IntoBpIterator: Copy {
    const BASES_PER_BYTE: usize;

    fn bp(&self) -> usize;

    fn sub_slice(&self, idx: usize, len: usize) -> Self;

    fn as_ptr(&self) -> *const u8;

    /// A simple iterator over characters.
    /// Returns u8 values in [0, 4).
    fn iter_bp(self) -> impl ExactSizeIterator<Item = u8>;

    /// An iterator that splits the input into 8 chunks and streams over them in parallel.
    /// Returns a separate `tail` iterator over the remaining characters.
    /// k is the k-mer size being iterated. When k>1, consecutive chunk overlap by k-1 characters.
    fn par_iter_bp(self, k: usize) -> (impl ExactSizeIterator<Item = S>, Self);
}

impl<'s> IntoBpIterator for &'s [u8] {
    const BASES_PER_BYTE: usize = 1;

    fn bp(&self) -> usize {
        self.len()
    }

    fn sub_slice(&self, idx: usize, len: usize) -> Self {
        &self[idx..idx + len]
    }

    fn as_ptr(&self) -> *const u8 {
        (self as &[u8]).as_ptr()
    }

    fn iter_bp(self) -> impl ExactSizeIterator<Item = u8> {
        self.iter().copied()
    }

    fn par_iter_bp(self, k: usize) -> (impl ExactSizeIterator<Item = S>, Self) {
        let num_kmers = self.len().saturating_sub(k - 1);
        let n = num_kmers / L;

        let offsets_lanes_0_4: Simd<*const u8, 4> =
            from_fn(|l| unsafe { self.as_ptr().add(l * n) }).into();
        let offsets_lanes_4_8: Simd<*const u8, 4> =
            from_fn(|l| unsafe { self.as_ptr().add((4 + l) * n) }).into();
        let mut upcoming_1 = S::splat(0);
        let mut upcoming_2 = S::splat(0);

        let it = (0..if num_kmers == 0 { 0 } else { n + k - 1 }).map(move |i| {
            if i % 4 == 0 {
                if i % 8 == 0 {
                    // Read a u64 containing the next 8 characters.
                    let idx_0_4 = offsets_lanes_0_4.wrapping_add(Simd::splat(i));
                    let idx_4_8 = offsets_lanes_4_8.wrapping_add(Simd::splat(i));
                    let u64_0_4: Simd<u32, 8> =
                        unsafe { transmute(Simd::<u64, 4>::gather_ptr(transmute(idx_0_4))) };
                    let u64_4_8: Simd<u32, 8> =
                        unsafe { transmute(Simd::<u64, 4>::gather_ptr(transmute(idx_4_8))) };
                    // Split into two vecs containing a u32 of 4 characters each.
                    (upcoming_1, upcoming_2) = u64_0_4.deinterleave(u64_4_8);
                } else {
                    // Move on to the next u32 containing 4 buffered characters.
                    upcoming_1 = upcoming_2;
                }
            }
            // Extract the last 2 bits of each character.
            let chars = upcoming_1 & S::splat(0x03);
            // Shift remaining characters to the right.
            upcoming_1 >>= S::splat(8);
            chars.into()
        });

        (it, &self[L * n..])
    }
}

#[derive(Copy, Clone)]
pub struct Packed<'s> {
    pub seq: &'s [u8],
    pub len_in_bp: usize,
}

impl<'s> IntoBpIterator for Packed<'s> {
    const BASES_PER_BYTE: usize = 4;

    fn bp(&self) -> usize {
        self.len_in_bp
    }

    fn sub_slice(&self, idx: usize, len: usize) -> Self {
        assert!(idx + len <= self.len_in_bp);
        assert!(idx % Self::BASES_PER_BYTE == 0);
        let seq = &self.seq[idx / 4..];
        Packed {
            seq,
            len_in_bp: len,
        }
    }

    fn as_ptr(&self) -> *const u8 {
        self.seq.as_ptr()
    }

    fn iter_bp(self) -> impl ExactSizeIterator<Item = u8> {
        assert!(self.len_in_bp <= self.seq.len() * 4);

        // read u64 at a time?
        let mut byte = 0;
        (0..self.len_in_bp).map(move |i| {
            if i % 4 == 0 {
                byte = self.seq[i / 4];
            }
            // Shift byte instead of i?
            (byte >> (2 * (i % 4))) & 0b11
        })
    }

    fn par_iter_bp(self, k: usize) -> (impl ExactSizeIterator<Item = S>, Self) {
        #[cfg(target_endian = "big")]
        assert!(false, "Big endian architectures are not yet supported.");

        let num_kmers = self.len_in_bp.saturating_sub(k - 1);
        let n = (num_kmers / L).prev_multiple_of(&4);
        let bytes_per_chunk = n / 4;

        let offsets_lanes_0_4: Simd<*const u8, 4> =
            from_fn(|l| unsafe { self.seq.as_ptr().add(l * bytes_per_chunk) }).into();
        let offsets_lanes_4_8: Simd<*const u8, 4> =
            from_fn(|l| unsafe { self.seq.as_ptr().add((4 + l) * bytes_per_chunk) }).into();
        let mut upcoming_1 = S::splat(0);
        let mut upcoming_2 = S::splat(0);

        let it = (0..if num_kmers == 0 { 0 } else { n + k - 1 }).map(move |i| {
            if i % 16 == 0 {
                if i % 32 == 0 {
                    // Read a u64 containing the next 8 characters.
                    let idx_0_4 = offsets_lanes_0_4.wrapping_add(Simd::splat(i / 4));
                    let idx_4_8 = offsets_lanes_4_8.wrapping_add(Simd::splat(i / 4));
                    let u64_0_4: Simd<u32, 8> =
                        unsafe { transmute(Simd::<u64, 4>::gather_ptr(transmute(idx_0_4))) };
                    let u64_4_8: Simd<u32, 8> =
                        unsafe { transmute(Simd::<u64, 4>::gather_ptr(transmute(idx_4_8))) };
                    // Split into two vecs containing a u32 of 4 characters each.
                    (upcoming_1, upcoming_2) = u64_0_4.deinterleave(u64_4_8);
                } else {
                    // Move on to the next u32 containing 4 buffered characters.
                    upcoming_1 = upcoming_2;
                }
            }
            // Extract the last 2 bits of each character.
            let chars = upcoming_1 & S::splat(0x03);
            // Shift remaining characters to the right.
            upcoming_1 >>= S::splat(2);
            chars.into()
        });

        (
            it,
            Packed {
                seq: &self.seq[L * bytes_per_chunk..],
                len_in_bp: self.len_in_bp - L * n,
            },
        )
    }
}
