use std::mem::transmute;

use packed_seq::Seq;
use wide::{i32x8, CmpGt};

use crate::Captures;

/// An iterator over windows that returns for each whether it's canonical or not.
/// Canonical windows have >half odd characters.
/// Window length l=k+w-1 must be odd for this to never tie.
pub fn canonical_scalar_it<'s>(
    seq: impl Seq<'s>,
    k: usize,
    w: usize,
) -> impl ExactSizeIterator<Item = bool> + Captures<&'s ()> {
    let l = k + w - 1;
    assert!(
        l % 2 == 1,
        "Window length must be odd to guarantee canonicality"
    );

    let mut add = seq.iter_bp();
    let remove = seq.iter_bp();

    // Cnt of odd characters, offset by -l/2 so >0 is canonical and <0 is not.
    let mut cnt = -(l as isize) / 2;

    add.by_ref().take(l - 1).for_each(|a| {
        cnt += a as isize & 1;
    });
    add.zip(remove).map(move |(a, r)| {
        cnt += a as isize & 1;
        cnt -= r as isize & 1;
        cnt > 0
    })
}

/// An iterator over windows that returns for each whether it's canonical or not.
/// Canonical windows have >half odd characters.
/// Window length l=k+w-1 must be odd for this to never tie.
///
/// Split the kmers of the sequence into 8 chunks of equal length ~len/8.
/// Then compute of each of them in parallel using SIMD,
/// and return the remaining few using the second iterator.
pub fn canonical_par_it<'s>(
    seq: impl Seq<'s>,
    k: usize,
    w: usize,
) -> (
    impl ExactSizeIterator<Item = i32x8> + Captures<&'s ()>,
    impl ExactSizeIterator<Item = bool> + Captures<&'s ()>,
) {
    let l = k + w - 1;
    assert!(
        l % 2 == 1,
        "Window length must be odd to guarantee canonicality"
    );

    // Cnt of odd characters, offset by -l/2 so >0 is canonical and <0 is not.
    let mut cnt = i32x8::splat(-(l as i32) / 2);

    let (mut add_remove, tail) = seq.par_iter_bp_delayed(k + w - 1, l - 1);

    let one = i32x8::splat(1);

    add_remove.by_ref().take(l - 1).for_each(|(a, _r)| {
        cnt += unsafe { transmute::<_, i32x8>(a) } & one;
    });

    let it = add_remove.map(
        #[inline(always)]
        move |(a, r)| {
            cnt += unsafe { transmute::<_, i32x8>(a) } & one;
            cnt -= unsafe { transmute::<_, i32x8>(r) } & one;
            cnt.cmp_gt(i32x8::splat(0))
        },
    );

    let tail = canonical_scalar_it(tail, k, w);

    (it, tail)
}