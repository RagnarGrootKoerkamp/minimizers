//! This code is adapted from the nthash crate:
use super::*;

#[derive(Clone, Copy, Debug)]
pub struct NtHash;
impl Hasher for NtHash {
    type Out = u64;
    #[inline(always)]
    fn hash(&self, t: &[u8]) -> u64 {
        ::nthash::ntf64(t, 0, t.len())
    }
    #[inline(always)]
    fn hash_kmers(&self, k: usize, t: &[u8]) -> impl Iterator<Item = Self::Out> {
        NtHashForwardIterator::new(t, k).unwrap()
    }
}

pub(crate) const MAXIMUM_K_SIZE: usize = u32::max_value() as usize;

const H_LOOKUP: [u64; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x3c8b_fbb3_95c6_0474;
    lookup[b'C' as usize] = 0x3193_c185_62a0_2b4c;
    lookup[b'G' as usize] = 0x2032_3ed0_8257_2324;
    lookup[b'T' as usize] = 0x2955_49f5_4be2_4456;
    lookup[b'N' as usize] = 0;
    lookup
};

#[inline(always)]
fn h(c: u8) -> u64 {
    H_LOOKUP[c as usize]
}

#[derive(Debug)]
pub struct NtHashForwardIterator<'a> {
    seq: &'a [u8],
    k: usize,
    fh: u64,
    current_idx: usize,
    max_idx: usize,
}

impl<'a> NtHashForwardIterator<'a> {
    /// Creates a new NtHashForwardIterator with internal state properly initialized.
    #[inline(always)]
    pub fn new(seq: &'a [u8], k: usize) -> Option<NtHashForwardIterator<'a>> {
        if k > seq.len() {
            return None;
        }
        if k > MAXIMUM_K_SIZE {
            return None;
        }

        let mut fh = 0;
        for (i, v) in seq[0..k].iter().enumerate() {
            fh ^= h(*v).rotate_left((k - i - 1) as u32);
        }

        Some(NtHashForwardIterator { seq, k, fh, current_idx: 0, max_idx: seq.len() - k + 1 })
    }
}

impl<'a> Iterator for NtHashForwardIterator<'a> {
    type Item = u64;

    #[inline(always)]
    fn next(&mut self) -> Option<u64> {
        if self.current_idx == self.max_idx {
            return None;
        };

        if self.current_idx != 0 {
            let i = self.current_idx - 1;
            let seqi = unsafe { *self.seq.get_unchecked(i) };
            let seqk = unsafe { *self.seq.get_unchecked(i + self.k) };

            self.fh = self.fh.rotate_left(1) ^ h(seqi).rotate_left(self.k as u32) ^ h(seqk);
        }

        self.current_idx += 1;
        Some(self.fh)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.max_idx, Some(self.max_idx))
    }
}

impl<'a> ExactSizeIterator for NtHashForwardIterator<'a> {}
