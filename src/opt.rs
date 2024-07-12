use std::{
    array::from_fn,
    simd::{cmp::SimdOrd, LaneCount, Simd, SimdElement, SupportedLaneCount},
};

use super::*;

/// Rehashing, based on Daniel Liu's implementation at
/// https://gist.github.com/Daniel-Liu-c0deb0t/7078ebca04569068f15507aa856be6e8
/// Idea: skip the priority queue and simply only track the smallest minimizer in the window.
/// When the smallest one falls out, simply re-scan the entire window.
pub struct MinimizerRehash<O: Order> {
    pub k: usize,
    pub w: usize,
    pub o: O,
}

impl<O: Order> SamplingScheme for MinimizerRehash<O> {
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        let mut min = usize::MAX;
        let mut min_idx = 0;
        let mut kmers = text.windows(self.k).enumerate();
        // Process the first w-1 kmers.
        for (j, kmer) in kmers.by_ref().take(self.w - 1) {
            let h = self.o.key(kmer);
            if h < min {
                min = h;
                min_idx = j;
            }
        }
        // i: absolute position of lmer
        // j: absolute position of kmer
        kmers.enumerate().map(move |(i, (j, kmer))| {
            if min_idx < i {
                min = usize::MAX;

                // Re-scan the window.
                for (j, kmer) in text.windows(self.k).enumerate().skip(i).take(self.w) {
                    let h = self.o.key(kmer);
                    if h < min {
                        min = h;
                        min_idx = j;
                    }
                }
            } else {
                let h = self.o.key(kmer);
                // Shift in the new one.
                if h < min {
                    min = h;
                    min_idx = j;
                }
            }
            min_idx
        })
    }
}

pub struct MinimizerRescan<O: Order> {
    pub k: usize,
    pub w: usize,
    pub o: O,
}

impl<O: Order> SamplingScheme for MinimizerRescan<O> {
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        let mut min = usize::MAX;
        let mut min_idx = 0;
        let mut kmers = text.windows(self.k).enumerate();

        // Rolling window of last w hashes.
        let mut hashes = vec![0; self.w];

        // Process the first w-1 kmers.
        for (j, kmer) in kmers.by_ref().take(self.w - 1) {
            let h = self.o.key(kmer);
            hashes[j] = h;
            if h < min {
                min = h;
                min_idx = j;
            }
        }
        let mut hash_idx = self.w - 1;

        // i: absolute position of lmer
        // j: absolute position of kmer
        kmers.enumerate().map(move |(i, (j, kmer))| {
            let h = self.o.key(kmer);
            hashes[hash_idx] = h;
            let cur_hash_idx = hash_idx;
            hash_idx += 1;
            if hash_idx == self.w {
                hash_idx = 0;
            }

            // w = 5
            //           j
            //  10  11  12 | 8   9
            //   0   1   2   3   4
            // |   |   | c | h |   |
            //               ^ new hash_idx
            //           ^ just inserted, cur_hash_idx

            if min_idx < i {
                // Re-scan the window.
                let p1 = hashes[cur_hash_idx + 1..]
                    .iter()
                    .position_min()
                    .map_or((usize::MAX, usize::MAX), |pos| {
                        (hashes[cur_hash_idx + 1 + pos], j + pos + 1 - self.w)
                    });
                let p2 = hashes[..=cur_hash_idx]
                    .iter()
                    .position_min()
                    .map_or((usize::MAX, usize::MAX), |pos| {
                        (hashes[pos], j - (cur_hash_idx - pos))
                    });
                (min, min_idx) = std::cmp::min(p2, p1);
            } else {
                // Shift in the new one.
                if h < min {
                    min = h;
                    min_idx = j;
                }
            }

            min_idx
        })
    }
}

pub trait HashPos: Clone + Copy {
    const MAX: Self;
    fn new(hash: u64, pos: usize) -> Self;
    fn pos(&self) -> usize;
    fn min(&self, other: &Self) -> Self;
    fn min_with(&mut self, other: &Self)
    where
        Self: Sized,
    {
        *self = self.min(other);
    }
    fn update_with(&mut self, hash: u64, pos: usize);
}

pub type Tuple = (u64, usize);
pub type Pack = u64;
pub type SignedPack = i64;

impl HashPos for (u64, usize) {
    const MAX: Self = (u64::MAX, usize::MAX);
    #[inline(always)]
    fn new(hash: u64, pos: usize) -> Self {
        (hash, pos)
    }

    #[inline(always)]
    fn pos(&self) -> usize {
        self.1
    }

    #[inline(always)]
    fn min(&self, other: &Self) -> Self {
        if *other < *self {
            *other
        } else {
            *self
        }
    }

    #[inline(always)]
    fn update_with(&mut self, hash: u64, pos: usize) {
        if hash < self.0 {
            self.0 = hash;
            self.1 = pos;
        }
    }
}

impl HashPos for u64 {
    const MAX: Self = u64::MAX;
    #[inline(always)]
    fn new(hash: u64, pos: usize) -> Self {
        (hash << 32) | pos as u64
    }

    #[inline(always)]
    fn pos(&self) -> usize {
        *self as u32 as usize
    }

    #[inline(always)]
    fn min(&self, other: &Self) -> Self {
        if *other < *self {
            *other
        } else {
            *self
        }
    }

    #[inline(always)]
    fn update_with(&mut self, hash: u64, pos: usize) {
        if (hash >> 32) < (*self >> 32) {
            *self = Self::new(hash, pos);
        }
    }
}

impl HashPos for i64 {
    const MAX: Self = i64::MAX;
    #[inline(always)]
    fn new(hash: u64, pos: usize) -> Self {
        ((hash << 32) | pos as u64) as i64
    }

    #[inline(always)]
    fn pos(&self) -> usize {
        *self as u32 as usize
    }

    #[inline(always)]
    fn min(&self, other: &Self) -> Self {
        if *other < *self {
            *other
        } else {
            *self
        }
    }

    #[inline(always)]
    fn update_with(&mut self, _hash: u64, _pos: usize) {
        unimplemented!();
    }
}

impl HashPos for u32 {
    const MAX: Self = u32::MAX;
    #[inline(always)]
    fn new(hash: u64, pos: usize) -> Self {
        ((hash << 16) | (pos as u64) & 0xffff) as u32
    }

    #[inline(always)]
    fn pos(&self) -> usize {
        (*self & 0xffff) as usize
    }

    #[inline(always)]
    fn min(&self, other: &Self) -> Self {
        if *other < *self {
            *other
        } else {
            *self
        }
    }

    #[inline(always)]
    fn update_with(&mut self, _hash: u64, _pos: usize) {
        unimplemented!();
    }
}

pub struct MinimizerRescanNt<T: HashPos = Tuple> {
    pub k: usize,
    pub w: usize,
    pub _t: std::marker::PhantomData<T>,
}

impl<T: HashPos> MinimizerRescanNt<T> {
    pub fn new(k: usize, w: usize) -> Self {
        Self {
            k,
            w,
            _t: std::marker::PhantomData,
        }
    }
}

impl<T: HashPos> SamplingScheme for MinimizerRescanNt<T> {
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        let mut min = T::MAX;
        let kmer_hashes = nthash::NtHashForwardIterator::new(text, self.k).unwrap();
        let mut kmers = kmer_hashes.enumerate();

        // Rolling window of last w hashes.
        let mut hashes = vec![0; self.w];

        // Process the first w-1 kmers.
        for (j, h) in kmers.by_ref().take(self.w - 1) {
            hashes[j] = h;
            let h = T::new(h, j);
            min.min_with(&h);
        }
        let mut hash_idx = self.w - 1;

        // i: absolute position of lmer
        // j: absolute position of kmer
        kmers.enumerate().map(
            #[inline(always)]
            move |(i, (j, h))| {
                unsafe { *hashes.get_unchecked_mut(hash_idx) = h };
                let cur_hash_idx = hash_idx;
                hash_idx += 1;
                if hash_idx == self.w {
                    hash_idx = 0;
                }

                if min.pos() >= i {
                    // Shift in the new one.
                    min.update_with(h, j);
                } else {
                    // The old minimum is dropped; re-scan the window.
                    let p1 =
                        hashes[cur_hash_idx + 1..]
                            .iter()
                            .position_min()
                            .map_or(T::MAX, |pos| {
                                T::new(
                                    unsafe { *hashes.get_unchecked(cur_hash_idx + 1 + pos) },
                                    j + pos + 1 - self.w,
                                )
                            });
                    let p2 = hashes[..=cur_hash_idx]
                        .iter()
                        .position_min()
                        .map_or(T::MAX, |pos| {
                            T::new(
                                unsafe { *hashes.get_unchecked(pos) },
                                j - (cur_hash_idx - pos),
                            )
                        });
                    min = p1.min(&p2);
                }

                min.pos()
            },
        )
    }
}

/// Idea based on https://codeforces.com/blog/entry/71687:
/// Every w positions, compute suffix minimimums for the w suffixes ending there.
pub struct MinimizerStacks<T: HashPos = Tuple> {
    pub k: usize,
    pub w: usize,
    pub _t: std::marker::PhantomData<T>,
}

impl<T: HashPos> MinimizerStacks<T> {
    pub fn new(k: usize, w: usize) -> Self {
        Self {
            k,
            w,
            _t: std::marker::PhantomData,
        }
    }
}

impl<T: HashPos> SamplingScheme for MinimizerStacks<T> {
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        let kmer_hashes = nthash::NtHashForwardIterator::new(text, self.k).unwrap();
        let mut kmers = kmer_hashes.enumerate();

        // Rolling window of last w hashes.
        let mut hashes = vec![T::MAX; self.w];

        // Process the first w-1 kmers.
        for (j, h) in kmers.by_ref().take(self.w - 1) {
            hashes[j] = T::new(h, j);
        }
        let mut hash_idx = self.w - 1;

        let mut rmin = T::MAX;

        // i: absolute position of lmer
        // j: absolute position of kmer
        kmers.map(
            #[inline(always)]
            move |(j, h)| {
                let h = T::new(h, j);
                unsafe { *hashes.get_unchecked_mut(hash_idx) = h };
                rmin.min_with(&h);
                hash_idx += 1;
                if hash_idx == self.w {
                    hash_idx = 0;

                    // Rolling suffix minima over the prefix.
                    for i in (0..self.w - 1).rev() {
                        unsafe {
                            let y = *hashes.get_unchecked(i + 1);
                            hashes.get_unchecked_mut(i).min_with(&y);
                        }
                    }

                    rmin = T::MAX;
                }
                unsafe { hashes.get_unchecked(hash_idx) }.min(&rmin).pos()
            },
        )
    }
}

/// Idea based on https://codeforces.com/blog/entry/71687:
/// Every w positions, compute suffix minimimums for the w suffixes ending there.
pub struct MinimizerStacksBuf<T: HashPos = Tuple> {
    pub k: usize,
    pub w: usize,
    pub _t: std::marker::PhantomData<T>,
}

impl<T: HashPos> MinimizerStacksBuf<T> {
    pub fn new(k: usize, w: usize) -> Self {
        Self {
            k,
            w,
            _t: std::marker::PhantomData,
        }
    }
}

impl<T: HashPos> SamplingScheme for MinimizerStacksBuf<T> {
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        let kmer_hashes = nthash::NtHashForwardIterator::new(text, self.k).unwrap();
        let mut kmers = kmer_hashes.enumerate();

        // Rolling window containing (pos hash, suffix hash, pos).
        let mut hashes = vec![(0, T::MAX); 3 * self.w];
        let mut i = self.w - 1;
        let mut read_l = 0 + i;
        let mut read_r = self.w + i;
        let mut update = 2 * self.w - 1 - i;
        let mut write = 2 * self.w + i;

        // Process the first w-1 kmers.
        for (j, h) in kmers.by_ref().take(self.w - 1) {
            hashes[2 * self.w + j] = (h, T::new(h, j));
        }

        // NOTE: The hash is in the upper 32 bits, and the position in the lower 32 bits.
        let mut rmin = T::MAX;
        let mut lmin = T::MAX;

        // i: absolute position of lmer
        // j: absolute position of kmer
        kmers.enumerate().map(
            #[inline(always)]
            move |(_i, (j, h))| {
                // write
                unsafe { *hashes.get_unchecked_mut(write) = (h, T::new(h, j)) };
                // update
                let entry = unsafe { &mut *hashes.get_unchecked_mut(update) };
                lmin = entry.1.min(&lmin);
                entry.1 = lmin;

                // read rmin
                let entry = unsafe { *hashes.get_unchecked(read_r) };
                rmin.min_with(&T::new(entry.0, j - 1 * self.w));

                i += 1;
                read_l += 1;
                read_r += 1;
                update -= 1;
                write += 1;
                if i == self.w {
                    i = 0;
                    (read_l, read_r, update, write) =
                        (read_r - self.w, write - self.w, write - 1, read_l - self.w);
                    lmin = T::MAX;
                    rmin = T::MAX;
                }

                // read lmin
                let lmin = unsafe { hashes.get_unchecked(read_l).1 };
                lmin.min(&rmin).pos()
            },
        )
    }
}

/////////////////////////////////////////////////////

pub struct MinimizerStacksSimd {
    pub k: usize,
    pub w: usize,
}

impl MinimizerStacksSimd {
    pub fn new(k: usize, w: usize) -> Self {
        assert!(w <= 256);
        Self { k, w }
    }
}

impl SamplingScheme for MinimizerStacksSimd {
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        type T = u32;
        const L: usize = 8;
        type S = Simd<T, L>;

        // split text in L chunks
        // TODO: Edge cases.
        let len = text.len() / L;
        // let text_chunks: [&[u8]; L] = text.chunks(len).collect_vec().try_into().unwrap();
        // let num_kmers = len - self.k + 1;

        let mut kmers = nthash::NtHashForwardIteratorSimd::new(text, len, self.k)
            .unwrap()
            .enumerate();

        // Rolling window of last w hashes.
        let mut hashes: Vec<S> = vec![S::splat(T::MAX); self.w];

        // Process the first w-1 kmers.
        kmers.by_ref().take(self.w - 1).for_each(|(j, h)| {
            hashes[j] = from_fn(|l| T::new(h[l] as u64, j)).into();
        });
        let mut hash_idx = self.w - 1;

        let mut rmin: S = S::splat(T::MAX);

        let pos_offset: Simd<u32, L> = from_fn(|l| (l * len) as u32).into();

        assert!(self.w > 1);

        // i: absolute position of lmer
        // j: absolute position of kmer
        let mut poss = vec![];
        poss.reserve(len);
        kmers.for_each(
            #[inline(always)]
            |(j, h)| {
                let h = from_fn(
                    #[inline(always)]
                    |l| T::new(h[l] as u64, j),
                )
                .into();

                // let h = h.as_array().map(|h| T::new(h as u64, j)).into();
                unsafe { *hashes.get_unchecked_mut(hash_idx) = h };
                rmin = rmin.simd_min(h);
                hash_idx += 1;
                if hash_idx == self.w {
                    hash_idx = 0;

                    // Rolling suffix minima over the prefix.
                    for i in (0..self.w - 1).rev() {
                        unsafe {
                            let y = *hashes.get_unchecked(i + 1);
                            let x = hashes.get_unchecked_mut(i);
                            *x = x.simd_min(y);
                        }
                    }

                    rmin = S::splat(T::MAX);
                }

                let ps: Simd<u32, L> = unsafe { hashes.get_unchecked(hash_idx) }
                    .simd_min(rmin)
                    .to_array()
                    .map(
                        #[inline(always)]
                        |l| l.pos() as u32,
                    )
                    .into();
                let ps = ps + pos_offset;
                poss.push(ps);
                // TODO: Gather 8-bit offsets into a u64 and scatter once every 8 iterations.
                // for l in 0..L {
                //     // poss[j + l * num_kmers] = ps[l];
                //     unsafe { *poss.get_unchecked_mut(j + l * num_kmers) = ps[l] };
                // }
            },
        );
        SimdIter::new(poss).map(
            #[inline(always)]
            |p| p as usize,
        )
    }
}

/// Iterates first over the first lane of all vector elements, then over the second lane, and so on.
struct SimdIter<T: SimdElement, const L: usize>
where
    LaneCount<L>: SupportedLaneCount,
{
    data: Vec<Simd<T, L>>,
    lane: usize,
    idx: usize,
}
impl<T: SimdElement, const L: usize> SimdIter<T, L>
where
    LaneCount<L>: SupportedLaneCount,
{
    #[inline(always)]
    fn new(data: Vec<Simd<T, L>>) -> Self {
        Self {
            data,
            lane: 0,
            idx: 0,
        }
    }
}

impl<T: SimdElement, const L: usize> Iterator for SimdIter<T, L>
where
    LaneCount<L>: SupportedLaneCount,
{
    type Item = T;
    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx == self.data.len() {
            cold();
            self.idx = 0;
            self.lane += 1;
            if self.lane == L {
                return None;
            }
        }
        let res = unsafe {
            *self
                .data
                .get_unchecked(self.idx)
                .as_array()
                .get_unchecked(self.lane)
        };
        self.idx += 1;
        Some(res)
    }
}

#[inline(never)]
#[cold]
fn cold() {}

// TODO: Precompute 4^2 lookup table -> not good with SIMD.
// TODO: Alternative hash: take xor of t=8-mers multiplied by constant C.
// TODO: https://en.algorithmica.org/hpc/algorithms/prefix/
