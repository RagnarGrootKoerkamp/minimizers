#![allow(dead_code)]
//! This file is mostly copied directly from the nthash crate, with modifications for higher performance.

//! ntHash is a hash function tuned for genomic data.
//! It performs best when calculating hash values for adjacent k-mers in
//! an input sequence, operating an order of magnitude faster than the best
//! performing alternatives in typical use cases.
//!
//! [Scientific article with more details](https://doi.org/10.1093/bioinformatics/btw397)
//!
//! [Original implementation in C++](https://github.com/bcgsc/ntHash/)
//!
//! This crate is based on ntHash [1.0.4](https://github.com/bcgsc/ntHash/releases/tag/v1.0.4).
//!

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
    let val = H_LOOKUP[c as usize];
    if val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", c as char)
    }
    val
}

/// Calculate the hash for a k-mer in the forward strand of a sequence.
///
/// This is a low level function, more useful for debugging than for direct use.
///
/// ```
///    use nthash::ntf64;
///    let fh = ntf64(b"TGCAG", 0, 5);
///    assert_eq!(fh, 0xbafa6728fc6dabf);
/// ```
pub fn ntf64(s: &[u8], i: usize, k: usize) -> u64 {
    let mut out = h(s[i + k - 1]);
    for (idx, v) in s.iter().skip(i).take(k - 1).enumerate() {
        out ^= h(*v).rotate_left((k - idx - 1) as u32);
    }
    out
}

/// An efficient iterator for calculating hashes for genomic sequences. This
/// returns the forward hashes, not the canonical hashes.
///
/// Since it implements the `Iterator` trait it also
/// exposes many other useful methods. In this example we use `collect` to
/// generate all hashes and put them in a `Vec<u64>`.
/// ```
///     # use nthash::Result;
///     use nthash::NtHashForwardIterator;
///
///     # fn main() -> Result<()> {
///     let seq = b"ACTGC";
///     let iter = NtHashForwardIterator::new(seq, 3)?;
///     let hashes: Vec<u64> = iter.collect();
///     assert_eq!(hashes, [0xb85d2431d9ba031e, 0xb4d7ab2f9f1306b8, 0xd4a29bf149877c5c]);
///     # Ok(())
///     # }
/// ```
/// or, in one line:
/// ```
///     # use nthash::Result;
///     use nthash::NtHashForwardIterator;
///
///     # fn main() -> Result<()> {
///     assert_eq!(NtHashForwardIterator::new(b"ACTGC", 3)?.collect::<Vec<u64>>(),
///                [0xb85d2431d9ba031e, 0xb4d7ab2f9f1306b8, 0xd4a29bf149877c5c]);
///     # Ok(())
///     # }
/// ```
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

        Some(NtHashForwardIterator {
            seq,
            k,
            fh,
            current_idx: 0,
            max_idx: seq.len() - k + 1,
        })
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
            let seqi = self.seq[i];
            let seqk = self.seq[i + self.k];

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
