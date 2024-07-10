use super::*;

/// Rehashing, based on Daniel Lius implementation at
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

pub struct MinimizerRescanNt {
    pub k: usize,
    pub w: usize,
}

impl SamplingScheme for MinimizerRescanNt {
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        let mut min = u64::MAX;
        let mut min_idx = 0;
        let kmer_hashes = nthash::NtHashForwardIterator::new(text, self.k).unwrap();
        let mut kmers = kmer_hashes.enumerate();

        // Rolling window of last w hashes.
        let mut hashes = vec![0; self.w];

        // Process the first w-1 kmers.
        for (j, h) in kmers.by_ref().take(self.w - 1) {
            hashes[j] = h;
            if h < min {
                min = h;
                min_idx = j;
            }
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

                if min_idx >= i {
                    // Shift in the new one.
                    if h < min {
                        min = h;
                        min_idx = j;
                    }
                } else {
                    // The old minimum is dropped; re-scan the window.
                    let p1 = hashes[cur_hash_idx + 1..].iter().position_min().map_or(
                        (u64::MAX, usize::MAX),
                        |pos| {
                            (
                                unsafe { *hashes.get_unchecked(cur_hash_idx + 1 + pos) },
                                j + pos + 1 - self.w,
                            )
                        },
                    );
                    let p2 = hashes[..=cur_hash_idx].iter().position_min().map_or(
                        (u64::MAX, usize::MAX),
                        |pos| {
                            (
                                unsafe { *hashes.get_unchecked(pos) },
                                j - (cur_hash_idx - pos),
                            )
                        },
                    );
                    (min, min_idx) = std::cmp::min_by_key(p1, p2, |x| x.0);
                }

                min_idx
            },
        )
    }
}

/// Idea based on https://codeforces.com/blog/entry/71687:
/// Every w positions, compute suffix minimimums for the w suffixes ending there.
pub struct MinimizerStacks {
    pub k: usize,
    pub w: usize,
}

impl SamplingScheme for MinimizerStacks {
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        let kmer_hashes = nthash::NtHashForwardIterator::new(text, self.k).unwrap();
        let mut kmers = kmer_hashes.enumerate();

        // Rolling window of last w hashes.
        let mut hashes = vec![(0, 0); self.w];

        // Process the first w-1 kmers.
        for (j, h) in kmers.by_ref().take(self.w - 1) {
            hashes[j] = (h, j);
        }
        let mut hash_idx = self.w - 1;

        let mut rmin = (u64::MAX, usize::MAX);

        // i: absolute position of lmer
        // j: absolute position of kmer
        kmers.enumerate().map(
            #[inline(always)]
            move |(_i, (j, h))| {
                unsafe { *hashes.get_unchecked_mut(hash_idx) = (h, j) };
                rmin = std::cmp::min_by_key(rmin, (h, j), |x| x.0);
                hash_idx += 1;
                if hash_idx == self.w {
                    hash_idx = 0;

                    // Rolling suffix minima over the prefix.
                    for i in (0..self.w - 1).rev() {
                        unsafe {
                            let y = *hashes.get_unchecked(i + 1);
                            let x = &mut *hashes.get_unchecked_mut(i);
                            *x = std::cmp::min_by_key(*x, y, |x| x.0);
                        }
                    }

                    rmin = (u64::MAX, usize::MAX);
                }
                std::cmp::min_by_key(unsafe { *hashes.get_unchecked(hash_idx) }, rmin, |x| x.0).1
            },
        )
    }
}

/// Idea based on https://codeforces.com/blog/entry/71687:
/// Every w positions, compute suffix minimimums for the w suffixes ending there.
pub struct MinimizerStacksBuf {
    pub k: usize,
    pub w: usize,
}

impl SamplingScheme for MinimizerStacksBuf {
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> impl MinimizerIt {
        let kmer_hashes = nthash::NtHashForwardIterator::new(text, self.k).unwrap();
        // .map(|x| x >> 58);
        let mut kmers = kmer_hashes.enumerate();

        // Rolling window containing (pos hash, suffix hash, pos).
        let mut hashes = vec![(0, 0, 0); 3 * self.w];
        let mut i = self.w - 1;
        let mut read_l = 0 + i;
        let mut read_r = self.w + i;
        let mut update = 2 * self.w - 1 - i;
        let mut write = 2 * self.w + i;

        // Process the first w-1 kmers.
        for (j, h) in kmers.by_ref().take(self.w - 1) {
            hashes[2 * self.w + j] = (h, h, j);
        }

        let mut rmin = (u64::MAX, usize::MAX);
        let mut lmin = (u64::MAX, usize::MAX);

        let cmp = |x: (u64, usize), y: (u64, usize)| std::cmp::min_by_key(x, y, |x| x.0);

        // i: absolute position of lmer
        // j: absolute position of kmer
        kmers.enumerate().map(
            #[inline(always)]
            move |(_i, (j, h))| {
                // write
                unsafe { *hashes.get_unchecked_mut(write) = (h, h, j) };
                // update
                let entry = unsafe { &mut *hashes.get_unchecked_mut(update) };
                lmin = cmp((entry.1, entry.2), lmin);
                (entry.1, entry.2) = lmin;

                // read rmin
                let entry = unsafe { *hashes.get_unchecked(read_r) };
                rmin = cmp(rmin, (entry.0, j - 1 * self.w));

                i += 1;
                read_l += 1;
                read_r += 1;
                update -= 1;
                write += 1;
                if i == self.w {
                    i = 0;
                    (read_l, read_r, update, write) =
                        (read_r - self.w, write - self.w, write - 1, read_l - self.w);
                    lmin.0 = u64::MAX;
                    rmin.0 = u64::MAX;
                }

                // read lmin
                let entry = unsafe { *hashes.get_unchecked(read_l) };
                let ans = cmp((entry.1, entry.2), rmin);

                ans.1
            },
        )
    }
}

// TODO: Buffered, doing the reverse loop w steps ahead.
// TODO: SIMD? Processing 4 chunks/sequences in parallel?
// TODO: Precompute 4^2 lookup table.
// TODO: Alternative hash: take xor of t=8-mers multiplied by constant C.
// TODO: https://en.algorithmica.org/hpc/algorithms/prefix/
