use super::*;

pub struct Split;

impl<V: Ord + Copy + Max> SlidingMin<V> for Split {
    fn sliding_min(&self, w: usize, it: impl Iterator<Item = V>) -> impl Iterator<Item = Elem<V>> {
        let mut prefix_min = Elem { val: V::MAX, pos: 0 };
        let mut ring_buf = RingBuf::new(w, prefix_min);

        it.enumerate()
            .map(move |(pos, val)| {
                let elem = Elem { val, pos };
                ring_buf.push(elem);
                prefix_min = prefix_min.min(elem);
                // After a chunk has been filled, compute suffix minima.
                if ring_buf.idx() == 0 {
                    let mut suffix_min = ring_buf[w - 1];
                    for elem in ring_buf[0..w - 1].iter_mut().rev() {
                        suffix_min = suffix_min.min(*elem);
                        *elem = suffix_min;
                    }
                    prefix_min = Elem { val: V::MAX, pos: 0 };
                }
                prefix_min.min(ring_buf[ring_buf.idx()])
            })
            .skip(w - 1)
    }
}
