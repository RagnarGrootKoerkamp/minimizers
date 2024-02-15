use std::collections::VecDeque;

/// A queue of (time, K) objects.
/// Both time and K values are always increasing, so that the smallest K value is always at the front.
/// When a K-value is pushed, strictly larger preceding values are removed.
/// Use `push_strict` to ensure strictly increasing elements.
pub struct MonotoneQueue<K: Ord> {
    q: VecDeque<(usize, K)>,
}

impl<K: Ord> MonotoneQueue<K>
where
    K: Copy,
{
    /// Initialize a new queue.
    pub fn new() -> Self {
        Self { q: VecDeque::new() }
    }

    /// Push element `k` with timestamp `t`.
    /// Strictly larger preceding `k` are removed, so that the queue is non-decreasing.
    pub fn push(&mut self, t: usize, k: K) {
        while let Some(&back) = self.q.back() {
            if back.1 > k {
                self.q.pop_back();
            } else {
                break;
            }
        }
        self.q.push_back((t, k));
    }

    /// Push element `k` with timestamp `t`.
    /// Larger or equal preceding `k` are removed, so that the queue is strictly increasing.
    pub fn push_strict(&mut self, t: usize, k: K) {
        while let Some(&back) = self.q.back() {
            if back.1 >= k {
                self.q.pop_back();
            } else {
                break;
            }
        }
        self.q.push_back((t, k));
    }

    /// Drop elements with timestamp <t.
    /// Then return the first element with timestamp >= t.
    pub fn pop(&mut self, t: usize) -> Option<&(usize, K)> {
        while let Some(&front) = self.q.front() {
            if front.0 < t {
                self.q.pop_front();
            } else {
                break;
            }
        }
        self.q.front()
    }
}
