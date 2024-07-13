use itertools::Itertools;

use super::*;

pub struct Rescan<Val: Ord> {
    w: usize,
    /// Position of next element.
    pos: usize,
    /// Index in `vals` of next element.
    idx: usize,
    /// Position of the smallest element in the last window.
    min_pos: usize,
    /// Value of the smallest element in the last window.
    min_val: Val,
    vals: Vec<Val>,
}

impl<Val: Ord + Copy + Max> SlidingMin<Val> for Rescan<Val> {
    #[inline(always)]
    fn new(w: usize) -> Self {
        assert!(w > 0);
        Self {
            w,
            pos: 0,
            idx: 0,
            min_pos: 0,
            min_val: Val::MAX,
            vals: vec![Val::MAX; w],
        }
    }

    #[inline(always)]
    fn push(&mut self, val: Val) -> Elem<Val> {
        unsafe {
            *self.vals.get_unchecked_mut(self.idx) = val;
            (self.min_val, self.min_pos) = if val < self.min_val {
                (val, self.pos)
            } else {
                (self.min_val, self.min_pos)
            };
            if self.pos - self.min_pos == self.w {
                // Find the position of the minimum, preferring older elements that
                // come *after* self.idx.
                let p1 = self
                    .vals
                    .get_unchecked(self.idx + 1..)
                    .iter()
                    .position_min();
                let p2 = self
                    .vals
                    .get_unchecked(..=self.idx)
                    .iter()
                    .position_min()
                    .unwrap();
                (self.min_val, self.min_pos) = if let Some(p1) = p1 {
                    let p1 = self.idx + 1 + p1;
                    if *self.vals.get_unchecked(p1) <= *self.vals.get_unchecked(p2) {
                        (
                            *self.vals.get_unchecked(p1),
                            self.pos - self.idx + p1 - self.w,
                        )
                    } else {
                        (*self.vals.get_unchecked(p2), self.pos - self.idx + p2)
                    }
                } else {
                    (*self.vals.get_unchecked(p2), self.pos - self.idx + p2)
                };
            }
        }

        self.pos += 1;
        self.idx += 1;
        if self.idx == self.w {
            self.idx = 0;
        }

        return Elem {
            pos: self.min_pos,
            val: self.min_val,
        };
    }
}

pub struct V4Rescan {
    pub w: usize,
    pub k: usize,
}

impl Minimizer for V4Rescan {
    fn minimizers(&self, text: &[u8]) -> Vec<usize> {
        let mut q = Rescan::new(self.w);
        let mut kmers = text.windows(self.k);
        // Inset the first w-1 k-mers, that do not yet form a full window.
        for kmer in kmers.by_ref().take(self.w - 1) {
            q.push(fxhash::hash(kmer));
        }
        kmers.map(|kmer| q.push(fxhash::hash(kmer)).pos).collect()
    }
}
