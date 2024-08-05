//! A fast implementation of ntHash.
//! Input:
//! - A 2-bit packed DNA sequence, only ACGT.
//! - It's length in bp.
//! - The k-mer length k.
//! Output:
//! - An iterator over all 32-hashes of all k-mers in the sequence.
//!
//! All functions assume that AVX2 is available.

// TODO: Write about 2bit encoded packed represenatation.

use super::packed::{IntoBpIterator, S};
use crate::par::packed::L;
use std::arch::x86_64::_mm256_permutevar_ps;

const HASHES: [u32; 4] = [
    0x3c8b_fbb3_95c6_0474u64 as u32,
    0x3193_c185_62a0_2b4cu64 as u32,
    0x2032_3ed0_8257_2324u64 as u32,
    0x2955_49f5_4be2_4456u64 as u32,
];

pub fn nthash32f_kmer(seq: impl IntoBpIterator) -> u32 {
    let mut h: u32 = 0;
    seq.iter_bp().for_each(|a| {
        h = h.rotate_left(1) ^ HASHES[a as usize];
    });
    h
}

pub fn nthash32f_scalar_it(
    seq: impl IntoBpIterator,
    k: usize,
) -> impl ExactSizeIterator<Item = u32> {
    assert!(k > 0);
    let mut h: u32 = 0;
    let mut add = seq.iter_bp();
    let remove = seq.iter_bp();
    add.by_ref().take(k - 1).for_each(|a| {
        h = h.rotate_left(1) ^ HASHES[a as usize];
    });
    add.zip(remove).map(move |(a, r)| {
        let h_out = h.rotate_left(1) ^ HASHES[a as usize];
        h = h_out ^ (HASHES[r as usize].rotate_left((k - 1) as u32));
        h_out
    })
}

// Split the kmers of the sequence into 8 chunks of equal length ~len/8.
// Then return the hashes of each of them, and the remaining few using the second iterator.
// The tail end has up to 31*8 = 248 elements.
// TODO: SMALL_K + reusing gathers?
// TODO: ARM variant
#[cfg(target_feature = "avx2")]
pub fn nthash32f_par_it(
    seq: impl IntoBpIterator,
    k: usize,
    w: usize,
) -> (
    impl ExactSizeIterator<Item = S>,
    impl ExactSizeIterator<Item = u32>,
) {
    assert!(k > 0);
    assert!(w > 0);
    // Each 128-bit half has a copy of the 4 32-bit hashes.
    let table: S = [0, 1, 2, 3, 0, 1, 2, 3].map(|c| HASHES[c as usize]).into();
    let table_rot: S = [0, 1, 2, 3, 0, 1, 2, 3]
        .map(|c| (HASHES[c as usize]).rotate_left(k as u32 - 1))
        .into();

    #[cfg(target_feature = "avx2")]
    {
        let table_lookup = |table: S, bps: S| -> S {
            use std::mem::transmute;
            unsafe { transmute(_mm256_permutevar_ps(transmute(table), transmute(bps))) }
        };

        let mut h = S::splat(0);
        let (mut add, tail) = seq.par_iter_bp(k + w - 1);
        let (remove, _tail) = seq.par_iter_bp(k + w - 1);

        add.by_ref().take(k - 1).for_each(|a| {
            h = ((h << 1) | (h >> 31)) ^ table_lookup(table, a);
        });

        let it = add.zip(remove).map(move |(a, r)| {
            let out = ((h << 1) | (h >> 31)) ^ table_lookup(table, a);
            h = out ^ table_lookup(table_rot, r);
            out
        });

        let tail = nthash32f_scalar_it(tail, k);

        (it, tail)
    }

    #[cfg(not(target_feature = "avx2"))]
    {
        static WARNED: LazyCell<()> = LazyCell::new(|| {
            eprintln!("nthash32f_par_it: AVX2 not available, falling back to scalar.");
        });

        *WARNED;

        let tail = seq.iter_bp(k + w - 1);
        let tail = nthash32f_scalar_it(tail, k);
        (empty(), tail)
    }
}

// This splits the input in ~2^15bp chunks, and computes the hashes of
// each chunk using 8 streams in parallel. Those results are buffered and a
// linear iterator over them is returned. Once it runs out, the next 2^15bp chunk is computed.
pub fn nthash32f_fast_it(seq: impl IntoBpIterator, k: usize) -> impl ExactSizeIterator<Item = u32> {
    NtHash32ForwardIterator::new(seq, k)
}

struct NtHash32ForwardIterator<BI: IntoBpIterator> {
    seq: BI,
    k: usize,

    chunk_size: usize,
    par_size: usize,
    num_chunks: usize,

    cache: Vec<S>,
    tail_cache: Vec<u32>,

    /// current chunk
    c: usize,
    /// current lane in chunk
    l: usize,
    /// current index
    i: usize,

    /// number of lanes in the current chunk.
    max_l: usize,
    /// length of the current lane
    max_i: usize,
}

impl<BI: IntoBpIterator> NtHash32ForwardIterator<BI> {
    fn new(seq: BI, k: usize) -> Self {
        assert!(k > 0);
        let len = seq.bp();
        let chunk_size = 1 << 13;
        assert!(chunk_size % L == 0);
        let par_size = chunk_size / L;

        let num_chunks = len.saturating_sub(k - 1).div_ceil(chunk_size);

        let mut cache = vec![S::splat(0); par_size];
        let mut tail_cache = vec![];
        if num_chunks == 1 {
            buffer_chunk(seq, k, &mut cache, &mut tail_cache);
        }

        Self {
            seq,
            k,
            cache,
            tail_cache,

            chunk_size,
            par_size,
            num_chunks,

            c: usize::MAX,
            i: par_size - 1,
            l: L - 1,

            max_i: par_size,
            max_l: L,
        }
    }
}

impl<BI: IntoBpIterator> Iterator for NtHash32ForwardIterator<BI> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        self.i += 1;
        if self.i == self.max_i {
            self.i = 0;
            self.l += 1;
            if self.l == self.max_l {
                self.l = 0;
                self.c = self.c.wrapping_add(1);
                if self.c < self.num_chunks {
                    let offset = self.c * self.chunk_size;
                    buffer_chunk(
                        self.seq.sub_slice(
                            offset,
                            (self.chunk_size + self.k - 1).min(self.seq.bp() - offset),
                        ),
                        self.k,
                        &mut self.cache,
                        &mut self.tail_cache,
                    );
                    if self.c + 1 < self.num_chunks {
                        assert_eq!(self.cache.len(), self.par_size);
                        assert!(self.tail_cache.is_empty());
                    } else {
                        assert!(self.cache.len() <= self.par_size);
                        self.max_i = self.cache.len();
                        if self.max_i == 0 {
                            self.c += 1;
                        }
                    }
                }
                if self.c == self.num_chunks {
                    // Run one extra iteration on the tail.
                    self.max_l = 1;
                    self.max_i = self.tail_cache.len();
                    assert!(self.tail_cache.len() <= self.par_size);
                    self.cache.resize(self.max_i, S::splat(0));
                    for i in 0..self.max_i {
                        self.cache[i][0] = self.tail_cache[i];
                    }
                    if self.max_i == 0 {
                        self.c += 1;
                    }
                }
                if self.c == self.num_chunks + 1 {
                    return None;
                }
            }
        }

        Some(self.cache[self.i][self.l])
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = if self.c >= self.num_chunks + 1 {
            0
        } else if self.c == self.num_chunks {
            self.max_i - self.i
        } else {
            let total_kmers = self.seq.bp() - self.k + 1;
            let processed_kmers = 1 + self.c * L * self.par_size + self.l * self.max_i + self.i;
            total_kmers - processed_kmers
        };
        (len, Some(len))
    }
}

impl<BI: IntoBpIterator> ExactSizeIterator for NtHash32ForwardIterator<BI> {}

fn buffer_chunk(seq: impl IntoBpIterator, k: usize, buf: &mut Vec<S>, tail_buf: &mut Vec<u32>) {
    let (mut par_it, mut tail) = nthash32f_par_it(seq, k, 1);
    buf.resize(par_it.size_hint().0, S::splat(0));
    for i in 0..buf.len() {
        unsafe {
            *buf.get_unchecked_mut(i) = par_it.next().unwrap();
        }
    }
    tail_buf.resize(tail.size_hint().0, 0);
    for i in 0..tail_buf.len() {
        unsafe {
            *tail_buf.get_unchecked_mut(i) = tail.next().unwrap();
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::par::packed::Packed;
    use rand::random;
    use std::{cell::LazyCell, iter::once};

    const BYTE_SEQ: LazyCell<Vec<u8>> =
        LazyCell::new(|| (0..1024).map(|_| random::<u8>() % 4).collect());

    const PACKED_SEQ: LazyCell<Vec<u8>> =
        LazyCell::new(|| (0..256).map(|_| random::<u8>()).collect());

    #[test]
    fn scalar_byte() {
        let seq = &**BYTE_SEQ;
        for k in [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65,
        ] {
            for len in (0..100).chain(once(1024)) {
                let seq = seq.sub_slice(0, len);
                let single = seq
                    .windows(k)
                    .map(|seq| nthash32f_kmer(seq))
                    .collect::<Vec<_>>();
                let scalar = nthash32f_scalar_it(seq, k).collect::<Vec<_>>();
                assert_eq!(single, scalar, "k={}, len={}", k, len);
            }
        }
    }

    #[test]
    fn parallel_byte() {
        let seq = &**BYTE_SEQ;
        for k in [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65,
        ] {
            for len in (0..100).chain(once(1024)) {
                let seq = seq.sub_slice(0, len);
                let scalar = nthash32f_scalar_it(seq, k).collect::<Vec<_>>();
                let parallel = nthash32f_fast_it(seq, k).collect::<Vec<_>>();
                assert_eq!(scalar, parallel, "k={}, len={}", k, len);
            }
        }
    }

    #[test]
    fn parallel_packed() {
        let seq = Packed {
            seq: &*PACKED_SEQ,
            len_in_bp: 1024,
        };
        for k in [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65,
        ] {
            for len in (0..100).chain(once(1024)) {
                let seq = seq.sub_slice(0, len);
                let scalar = nthash32f_scalar_it(seq, k).collect::<Vec<_>>();
                let parallel = nthash32f_fast_it(seq, k).collect::<Vec<_>>();
                assert_eq!(scalar, parallel, "k={}, len={}", k, len);
            }
        }
    }

    #[test]
    fn parallel_iter_byte() {
        let seq = &**BYTE_SEQ;
        for k in [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65,
        ] {
            for len in (0..100).chain(once(1024)) {
                let seq = seq.sub_slice(0, len);
                let scalar = nthash32f_scalar_it(seq, k).collect::<Vec<_>>();
                let (par_head, tail) = nthash32f_par_it(seq, k, 1);
                let par_head = par_head.collect::<Vec<_>>();
                let parallel_iter = (0..L)
                    .flat_map(|l| par_head.iter().map(move |x| x[l]))
                    .chain(tail)
                    .collect::<Vec<_>>();

                assert_eq!(scalar, parallel_iter, "k={}, len={}", k, len);
            }
        }
    }

    #[test]
    fn parallel_iter_packed() {
        let seq = Packed {
            seq: &*PACKED_SEQ,
            len_in_bp: 1024,
        };
        for k in [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65,
        ] {
            for len in (0..100).chain(once(1024)) {
                let seq = seq.sub_slice(0, len);
                let scalar = nthash32f_scalar_it(seq, k).collect::<Vec<_>>();
                let (par_head, tail) = nthash32f_par_it(seq, k, 1);
                let par_head = par_head.collect::<Vec<_>>();
                let parallel_iter = (0..L)
                    .flat_map(|l| par_head.iter().map(move |x| x[l]))
                    .chain(tail)
                    .collect::<Vec<_>>();
                assert_eq!(scalar, parallel_iter, "k={}, len={}", k, len);
            }
        }
    }
}
