use std::{array, mem::transmute};

use packed_seq::S;
use wide::u32x8;

use super::dedup::write_unique_with_old;

/// Collect a parallel stream into a single vector.
/// Works by taking 8 elements from each stream, and then transposing the SIMD-matrix before writing out the results.
#[inline(always)]
pub fn collect(
    par_head: impl ExactSizeIterator<Item = S>,
    tail: impl ExactSizeIterator<Item = u32>,
) -> Vec<u32> {
    let len = par_head.len();
    let mut v = vec![0; len * 8 + tail.len()];

    let mut m = [unsafe { transmute([0; 8]) }; 8];
    let mut i = 0;
    par_head.for_each(|x| {
        m[i % 8] = x;
        if i % 8 == 7 {
            let t = transpose(m);
            for j in 0..8 {
                unsafe {
                    *v.get_unchecked_mut(j * len + 8 * (i / 8)..j * len + 8 * (i / 8) + 8)
                        .split_first_chunk_mut::<8>()
                        .unwrap()
                        .0 = transmute(t[j]);
                }
            }
        }
        i += 1;
    });

    // Manually write the unfinished parts of length k=i%8.
    let t = transpose(m);
    let k = i % 8;
    for j in 0..8 {
        unsafe {
            v[j * len + 8 * (i / 8)..j * len + 8 * (i / 8) + k]
                .copy_from_slice(&transmute::<_, [u32; 8]>(t[j])[..k]);
        }
    }

    // Manually write the explicit tail.
    for (i, x) in tail.enumerate() {
        v[8 * len + i] = x;
    }

    v
}

/// A utility function for creating masks to use with Intel shuffle and
/// permute intrinsics.
///
/// Copied from the standard library, since it is unstable.
const fn _mm_shuffle(z: u32, y: u32, x: u32, w: u32) -> i32 {
    ((z << 6) | (y << 4) | (x << 2) | w) as i32
}

/// Transpose a matrix of 8 SIMD vectors.
/// https://stackoverflow.com/questions/25622745/transpose-an-8x8-float-using-avx-avx2
///
/// TODO: Investigate other transpose functions mentioned there?
#[inline(always)]
fn transpose(m: [S; 8]) -> [S; 8] {
    unsafe {
        use std::arch::x86_64::*;
        let m: [__m256; 8] = transmute(m);
        let t0 = _mm256_unpacklo_ps(m[0], m[1]);
        let t1 = _mm256_unpackhi_ps(m[0], m[1]);
        let t2 = _mm256_unpacklo_ps(m[2], m[3]);
        let t3 = _mm256_unpackhi_ps(m[2], m[3]);
        let t4 = _mm256_unpacklo_ps(m[4], m[5]);
        let t5 = _mm256_unpackhi_ps(m[4], m[5]);
        let t6 = _mm256_unpacklo_ps(m[6], m[7]);
        let t7 = _mm256_unpackhi_ps(m[6], m[7]);
        let tt0 = _mm256_shuffle_ps(t0, t2, _mm_shuffle(1, 0, 1, 0));
        let tt1 = _mm256_shuffle_ps(t0, t2, _mm_shuffle(3, 2, 3, 2));
        let tt2 = _mm256_shuffle_ps(t1, t3, _mm_shuffle(1, 0, 1, 0));
        let tt3 = _mm256_shuffle_ps(t1, t3, _mm_shuffle(3, 2, 3, 2));
        let tt4 = _mm256_shuffle_ps(t4, t6, _mm_shuffle(1, 0, 1, 0));
        let tt5 = _mm256_shuffle_ps(t4, t6, _mm_shuffle(3, 2, 3, 2));
        let tt6 = _mm256_shuffle_ps(t5, t7, _mm_shuffle(1, 0, 1, 0));
        let tt7 = _mm256_shuffle_ps(t5, t7, _mm_shuffle(3, 2, 3, 2));
        let mut m: [__m256; 8] = [transmute([0; 8]); 8];
        m[0] = _mm256_permute2f128_ps(tt0, tt4, 0x20);
        m[1] = _mm256_permute2f128_ps(tt1, tt5, 0x20);
        m[2] = _mm256_permute2f128_ps(tt2, tt6, 0x20);
        m[3] = _mm256_permute2f128_ps(tt3, tt7, 0x20);
        m[4] = _mm256_permute2f128_ps(tt0, tt4, 0x31);
        m[5] = _mm256_permute2f128_ps(tt1, tt5, 0x31);
        m[6] = _mm256_permute2f128_ps(tt2, tt6, 0x31);
        m[7] = _mm256_permute2f128_ps(tt3, tt7, 0x31);
        transmute(m)
    }
}
