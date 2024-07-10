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
