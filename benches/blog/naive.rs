#![cfg_attr(any(), rustfmt::skip)]
use super::*;

pub struct Naive {
    pub w: usize,
    pub k: usize,
}

impl Minimizer for Naive {
    fn window_minimizers(&self, text: &[u8]) -> Vec<usize> {
        // Iterate over the windows of size l=w+k-1.
        text.windows(self.w + self.k - 1)
            .enumerate()
            // For each window, starting at pos j, find the lexicographically smallest k-mer.
            .map(|(j, window)| {
                j + window
                    .windows(self.k)
                    .enumerate()
                    .min_by_key(|(_idx, kmer)| fxhash::hash64(kmer))
                    .unwrap()
                    .0
            })
            .collect()
    }
}

pub struct BufferedSlidingMin;

impl<V: Copy + Max + Ord> SlidingMin<V> for BufferedSlidingMin {
    fn sliding_min(w: usize, it: impl Iterator<Item = V>) -> impl Iterator<Item = Elem<V>> {
        // Iterate over items and their positions.
        let it = it.enumerate().map(|(pos, val)| Elem { val, pos });
        // A ring buffer that holds the w last elements.
        let mut buf: Vec<_> = vec![Elem { val: V::MAX, pos: 0 }; w];
        // The index to write to. Always equals pos % w, but avoids the modulo.
        let mut idx = 0;
        it.map(move |elem| {
            buf[idx] = elem;
            idx += 1;
            if idx == w {
                idx = 0;
            }
            // Return the minimum element in the buffer.
            // Ties between values are broken by the position.
            *buf.iter().min().unwrap()
        })
        // The first w-1 elements do not complete a window.
        .skip(w - 1)
    }
}
