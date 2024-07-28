use super::*;

pub struct Rescan2;

impl<V: Ord + Copy + Max + std::fmt::Debug> SlidingMin<V> for Rescan2 {
    fn sliding_min(&self, w: usize, it: impl Iterator<Item = V>) -> impl Iterator<Item = Elem<V>> {
        let mut min = Elem { val: V::MAX, pos: 0 };
        let mut min2 = Elem { val: V::MAX, pos: 0 };
        let mut ring_buf = RingBuf::new(w, min);

        it.enumerate()
            .map(move |(pos, val)| {
                let elem = Elem { val, pos };
                ring_buf.push(elem);
                if elem < min2 {
                    if elem < min {
                        min2 = min;
                        min = elem
                    } else {
                        min2 = elem;
                    }
                }
                // If the minimum falls out of the window, rescan from min2.pos onward to find the new minimum.
                if pos - min.pos == w {
                    let range = min2.pos..=pos;
                    min = Elem { val: V::MAX, pos: 0 };
                    min2 = Elem { val: V::MAX, pos: 0 };
                    for pos in range.rev() {
                        let elem = ring_buf[pos % w];
                        if elem < min2 {
                            if elem < min {
                                min2 = min;
                                min = elem
                            } else {
                                min2 = elem;
                            }
                        }
                    }
                }
                min
            })
            .skip(w - 1)
    }
}
