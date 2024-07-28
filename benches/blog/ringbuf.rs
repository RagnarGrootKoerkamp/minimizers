pub struct RingBuf<V> {
    w: usize,
    idx: usize,
    data: Vec<V>,
}

impl<V: Clone> RingBuf<V> {
    pub fn idx(&self) -> usize {
        self.idx
    }

    pub fn new(w: usize, v: V) -> Self {
        let data = vec![v; w];
        RingBuf { w, idx: 0, data }
    }

    pub fn push(&mut self, v: V) {
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

    fn deref(&self) -> &[V] {
        &self.data
    }
}

impl<V> std::ops::DerefMut for RingBuf<V> {
    fn deref_mut(&mut self) -> &mut [V] {
        &mut self.data
    }
}
