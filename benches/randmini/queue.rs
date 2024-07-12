use std::collections::VecDeque;

use super::{Elem, Minimizer, SlidingMin};

pub struct MonotoneQueueA<Val: Ord> {
    w: usize,
    pos: usize,
    /// A queue of (pos, val) objects.
    /// Both pos and val values are always increasing, so that the smallest
    /// value is always at the front.
    q: VecDeque<Elem<Val>>,
}

impl<Val: Ord + Copy> SlidingMin<Val> for MonotoneQueueA<Val> {
    fn new(w: usize) -> Self {
        assert!(w > 0);
        Self {
            w,
            pos: 0,
            q: VecDeque::new(),
        }
    }

    fn push(&mut self, val: Val) -> Elem<Val> {
        // Strictly larger preceding `k` are removed, so that the queue remains
        // non-decreasing.
        while let Some(back) = self.q.back() {
            if back.val > val {
                self.q.pop_back();
            } else {
                break;
            }
        }
        self.q.push_back(Elem { pos: self.pos, val });
        let front = self.q.front().unwrap(); // Safe, because we just pushed.
        if self.pos - front.pos >= self.w {
            self.q.pop_front();
        }
        self.pos += 1;
        *self.q.front().unwrap() // Safe, because w > 0.
    }
}

pub struct V3AQueue {
    pub w: usize,
    pub k: usize,
}

impl Minimizer for V3AQueue {
    fn minimizers(&self, text: &[u8]) -> Vec<usize> {
        let mut q = MonotoneQueueA::new(self.w);
        let mut kmers = text.windows(self.k);
        // Inset the first w-1 k-mers, that do not yet form a full window.
        for kmer in kmers.by_ref().take(self.w - 1) {
            q.push(fxhash::hash(kmer));
        }
        kmers.map(|kmer| q.push(fxhash::hash(kmer)).pos).collect()
    }
}

pub struct MonotoneQueueB<Val: Ord> {
    w: usize,
    pos: usize,
    /// A queue of (pos, val) objects.
    /// Both pos and val values are always increasing, so that the smallest
    /// value is always at the front.
    q: VecDeque<Elem<Val>>,
}

impl<Val: Ord + Copy> SlidingMin<Val> for MonotoneQueueB<Val> {
    #[inline(always)]
    fn new(w: usize) -> Self {
        assert!(w > 0);
        Self {
            w,
            pos: 0,
            q: VecDeque::new(),
        }
    }

    #[inline(always)]
    fn push(&mut self, val: Val) -> Elem<Val> {
        // Strictly larger preceding `k` are removed, so that the queue remains
        // non-decreasing.
        while let Some(back) = self.q.back() {
            if back.val > val {
                self.q.pop_back();
            } else {
                break;
            }
        }
        self.q.push_back(Elem { pos: self.pos, val });
        let front = self.q.front().unwrap(); // Safe, because we just pushed.
        if self.pos - front.pos >= self.w {
            self.q.pop_front();
        }
        self.pos += 1;
        *self.q.front().unwrap() // Safe, because w > 0.
    }
}

pub struct V3BQueue {
    pub w: usize,
    pub k: usize,
}

impl Minimizer for V3BQueue {
    fn minimizers(&self, text: &[u8]) -> Vec<usize> {
        let mut q = MonotoneQueueB::new(self.w);
        let mut kmers = text.windows(self.k);
        // Inset the first w-1 k-mers, that do not yet form a full window.
        for kmer in kmers.by_ref().take(self.w - 1) {
            q.push(fxhash::hash(kmer));
        }
        kmers.map(|kmer| q.push(fxhash::hash(kmer)).pos).collect()
    }
}
