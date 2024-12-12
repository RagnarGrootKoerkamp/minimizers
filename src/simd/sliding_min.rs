use core::array::from_fn;
use packed_seq::S;
use std::hint::assert_unchecked;

/// A custom RingBuf implementation that has a fixed size `w` and wraps around.
struct RingBuf<V> {
    w: usize,
    idx: usize,
    data: Vec<V>,
}

impl<V: Clone> RingBuf<V> {
    #[inline(always)]
    fn new(w: usize, v: V) -> Self {
        assert!(w > 0);
        let data = vec![v; w];
        RingBuf { w, idx: 0, data }
    }

    /// Returns the next index to be written.
    #[inline(always)]
    fn idx(&self) -> usize {
        self.idx
    }

    #[inline(always)]
    fn push(&mut self, v: V) {
        *unsafe { self.data.get_unchecked_mut(self.idx) } = v;
        self.idx += 1;
        if self.idx == self.w {
            self.idx = 0;
        }
    }
}

/// A RingBuf can be used as a slice.
impl<V> std::ops::Deref for RingBuf<V> {
    type Target = [V];

    #[inline(always)]
    fn deref(&self) -> &[V] {
        &self.data
    }
}

/// A RingBuf can be used as a mutable slice.
impl<V> std::ops::DerefMut for RingBuf<V> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut [V] {
        &mut self.data
    }
}

#[inline(always)]
pub fn sliding_min_scalar_it(
    it: impl ExactSizeIterator<Item = u32>,
    w: usize,
) -> impl ExactSizeIterator<Item = u32> {
    assert!(w > 0);
    assert!(
        w < (1 << 15),
        "sliding_min is not tested for windows of length > 2^15."
    );
    assert!(
        it.size_hint().0 < (1 << 32),
        "sliding_min returns 32bit indices. Try splitting the input into 4GB chunks first."
    );
    let mut prefix_min = u32::MAX;
    let mut ring_buf = RingBuf::new(w, prefix_min);
    // We only compare the upper 16 bits of each hash.
    // Ties are broken automatically in favour of lower pos.
    let val_mask = 0xffff_0000;
    let pos_mask = 0x0000_ffff;
    let mut pos = 0;
    let max_pos = (1 << 16) - 1;
    let mut pos_offset = 0;

    let mut it = it.map(
        #[inline(always)]
        move |val| {
            // Make sure the position does not interfere with the hash value.
            if pos == max_pos {
                let delta = ((1 << 16) - 2 - w) as u32;
                pos -= delta;
                prefix_min -= delta;
                pos_offset += delta;
                for x in &mut *ring_buf {
                    *x -= delta;
                }
            }
            let elem = (val & val_mask) | pos;
            pos += 1;
            ring_buf.push(elem);
            prefix_min = prefix_min.min(elem);
            // After a chunk has been filled, compute suffix minima.
            if ring_buf.idx() == 0 {
                let mut suffix_min = ring_buf[w - 1];
                for i in (0..w - 1).rev() {
                    suffix_min = suffix_min.min(ring_buf[i]);
                    ring_buf[i] = suffix_min;
                }
                prefix_min = elem; // slightly faster than assigning S::splat(u32::MAX)
            }
            let suffix_min = unsafe { *ring_buf.get_unchecked(ring_buf.idx()) };
            (prefix_min.min(suffix_min) & pos_mask) + pos_offset
        },
    );
    // This optimizes better than it.skip(w-1).
    it.by_ref().take(w - 1).for_each(drop);
    it
}

#[inline(always)]
pub fn sliding_min_par_it(
    it: impl ExactSizeIterator<Item = S>,
    w: usize,
) -> impl ExactSizeIterator<Item = S> {
    assert!(w > 0);
    assert!(w < (1 << 15), "This method is not tested for large w.");
    assert!(it.size_hint().0 * 8 < (1 << 32));
    let mut prefix_min = S::splat(u32::MAX);
    let mut ring_buf = RingBuf::new(w, prefix_min);
    // We only compare the upper 16 bits of each hash.
    // Ties are broken automatically in favour of lower pos.
    let val_mask = S::splat(0xffff_0000);
    let pos_mask = S::splat(0x0000_ffff);
    let max_pos = S::splat((1 << 16) - 1);
    let mut pos = S::splat(0);
    let mut pos_offset: S =
        from_fn(|l| (l * (it.size_hint().0.saturating_sub(w - 1))) as u32).into();

    let mut it = it.map(
        #[inline(always)]
        move |val| {
            // Make sure the position does not interfere with the hash value.
            if pos == max_pos {
                // Slow case extracted to a function to have better inlining here.
                reset_positions_offsets(
                    w,
                    &mut pos,
                    &mut prefix_min,
                    &mut pos_offset,
                    &mut ring_buf,
                );
            }
            // slightly faster than assigning S::splat(u32::MAX)
            let elem = (val & val_mask) | pos;
            pos += S::splat(1);
            ring_buf.push(elem);
            prefix_min = prefix_min.min(elem);
            // After a chunk has been filled, compute suffix minima.
            if ring_buf.idx() == 0 {
                // Slow case extracted to a function to have better inlining here.
                suffix_minima(&mut ring_buf, w, &mut prefix_min, elem);
            }
            let suffix_min = unsafe { *ring_buf.get_unchecked(ring_buf.idx()) };
            (prefix_min.min(suffix_min) & pos_mask) + pos_offset
        },
    );
    // This optimizes better than it.skip(w-1).
    it.by_ref().take(w - 1).for_each(drop);
    it
}

fn suffix_minima(ring_buf: &mut RingBuf<S>, w: usize, prefix_min: &mut S, elem: S) {
    // Avoid some bounds checks when this function is not inlined.
    unsafe { assert_unchecked(ring_buf.len() == w) };
    unsafe { assert_unchecked(w > 0) };
    let mut suffix_min = ring_buf[w - 1];
    for i in (0..w - 1).rev() {
        suffix_min = suffix_min.min(ring_buf[i]);
        ring_buf[i] = suffix_min;
    }
    *prefix_min = elem;
}

fn reset_positions_offsets(
    w: usize,
    pos: &mut S,
    prefix_min: &mut S,
    pos_offset: &mut S,
    ring_buf: &mut RingBuf<S>,
) {
    let delta = S::splat((1 << 16) - 2 - w as u32);
    *pos -= delta;
    *prefix_min -= delta;
    *pos_offset += delta;
    for x in &mut **ring_buf {
        *x -= delta;
    }
}
