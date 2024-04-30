#![feature(exclusive_range_pattern, type_alias_impl_trait, iter_repeat_n)]

pub mod de_bruijn_seq;
pub mod ilp_scheme;
pub mod monotone_queue;
pub mod order;

use crate::{
    de_bruijn_seq::{cyclic_exact_density_string, exact_density_string},
    monotone_queue::MonotoneQueue,
};
use itertools::Itertools;
use order::*;
use rayon::iter::{ParallelBridge, ParallelIterator};
use std::{
    f32::consts::PI,
    iter::{repeat_n, zip},
};

pub trait Minimizers: Iterator<Item = usize> {}
impl<T: Iterator<Item = usize>> Minimizers for T {}

fn pack(kmer: &[u8], sigma: usize) -> usize {
    let mut v = 0;
    for c in kmer {
        v = sigma * v + *c as usize;
    }
    v
}

pub fn minimizer(s: &[u8], k: usize, o: &impl DirectedOrder) -> usize {
    assert!(k > 0);
    assert!(k <= s.len());
    s.windows(k)
        .enumerate()
        .min_by_key(|&(_i, kmer)| o.key(kmer))
        .unwrap()
        .0
}

pub fn random_minimizer(s: &[u8], k: usize) -> usize {
    minimizer(s, k, &order::RandomOrder)
}

/// TODO: Rolling hash using NtHash.
/// TODO: Take an iterator over u8 instead?
pub fn text_minimizers<'a>(
    text: &'a [u8],
    w: usize,
    k: usize,
    o: &'a impl DirectedOrder,
) -> impl Minimizers + 'a {
    let mut q = monotone_queue::MonotoneQueue::new();
    let mut kmers = text.windows(k).enumerate();
    for (j, kmer) in kmers.by_ref().take(w - 1) {
        q.push(j, o.key(kmer));
    }
    // i: position of lmer
    // j: position of kmer
    kmers.enumerate().map(move |(i, (j, kmer))| {
        q.push(j, o.key(kmer));
        q.pop(i).unwrap().0
    })
}

pub fn cyclic_minimizers<'a>(
    text: &'a [u8],
    len: usize,
    w: usize,
    k: usize,
    o: &'a impl DirectedOrder,
) -> usize {
    let mut q = monotone_queue::MonotoneQueue::new();
    let mut kmers = text.windows(k).take(len + w - 1).enumerate();
    for (j, kmer) in kmers.by_ref().take(w - 1) {
        q.push(j, o.key(kmer));
    }
    // i: position of lmer
    // j: position of kmer
    let mut poss = kmers
        .enumerate()
        .map(move |(i, (j, kmer))| {
            q.push(j, o.key(kmer));
            q.pop(i).unwrap().0 % len
        })
        .dedup()
        .collect_vec();
    poss.sort();
    poss.dedup();
    poss.len()
}

pub fn text_random_minimizers<'a>(text: &'a [u8], w: usize, k: usize) -> impl Minimizers + 'a {
    text_minimizers(text, w, k, &order::RandomOrder)
}

#[derive(Debug, Clone)]
pub struct LocalScheme {
    pub k: usize,
    pub w: usize,
    pub sigma: usize,
    /// The index in [w] to choose for each of sigma^(k+w-1) possible l-mers.
    pub map: Vec<u8>,
}

impl LocalScheme {
    pub fn get(&self, lmer: &[u8]) -> usize {
        debug_assert_eq!(self.map.len(), self.sigma.pow((self.k + self.w - 1) as u32));
        debug_assert_eq!(lmer.len(), self.k + self.w - 1);
        let mut v = 0;
        for c in lmer {
            assert!(*c < self.sigma as u8);
            v = self.sigma * v + *c as usize;
        }
        // eprintln!("lmer {lmer:?} sigma {} => v {v}", self.sigma);
        self.map[v] as usize
    }

    /// Check if the local scheme corresponds to a directed order on kmers.
    /// If so, return one such order.
    pub fn to_directed_order(&self) -> Option<ExplicitDirectedOrder> {
        let l = self.k + self.w - 1;
        let num_lmers = self.sigma.pow(l as u32);
        assert_eq!(self.map.len(), num_lmers);
        let num_kmers = self.sigma.pow(self.k as u32);
        // smaller[i] is a list of kmers j with i < j, indicating that i comes before j in the toposort.
        let mut smaller = vec![vec![]; num_kmers];
        // direction in which kmer i should be compared in case of equality.
        let mut directions = vec![None; num_kmers];

        for lmer in 0..num_lmers {
            let mut chars = vec![];
            {
                let mut lmer = lmer;
                for _ in 0..l {
                    chars.push((lmer % self.sigma) as u8);
                    lmer /= self.sigma;
                }
                chars.reverse();
            }

            let kmers = {
                let mut lmer = lmer;
                let mut kmers = (0..self.w)
                    .map(|_| {
                        let kmer = lmer % num_kmers;
                        lmer /= self.sigma;
                        kmer
                    })
                    .collect_vec();
                kmers.reverse();
                kmers
            };

            let chosen_idx = self.map[lmer] as usize;
            let chosen_kmer = kmers[chosen_idx];
            // eprintln!("lmer {lmer:>3}: {:?} => {chosen_idx}", chars);
            for i in 0..self.w {
                if i == chosen_idx {
                    continue;
                }
                let kmer = kmers[i];
                if kmer != chosen_kmer {
                    smaller[chosen_kmer].push(kmer);
                    // eprintln!(
                    //     "{:?} < {:?}",
                    //     &chars[chosen_idx..chosen_idx + self.k],
                    //     &chars[i..i + self.k]
                    // );
                } else {
                    let wanted_direction = if i < chosen_idx {
                        Some(Direction::Rightmost)
                    } else {
                        Some(Direction::Leftmost)
                    };
                    if directions[chosen_kmer].is_none() {
                        directions[chosen_kmer] = wanted_direction;
                    } else {
                        if directions[chosen_kmer] != wanted_direction {
                            eprintln!(
                                "Inconsistent direction for kmer {:?}",
                                &chars[i..i + self.k]
                            );
                            // return None;
                        }
                    }
                }
            }
        }

        for x in &mut smaller {
            x.sort();
            x.dedup();
        }

        // Toposort the 'smaller' graph.
        let mut order = vec![];
        {
            let mut visited = vec![0; num_kmers];
            fn dfs(
                i: usize,
                smaller: &Vec<Vec<usize>>,
                visited: &mut Vec<u8>,
                order: &mut Vec<usize>,
            ) -> Option<()> {
                if visited[i] == 1 {
                    eprintln!("Cycle detected");
                    return None;
                }
                if visited[i] == 2 {
                    return Some(());
                }
                visited[i] = 1;
                for &j in &smaller[i] {
                    dfs(j, smaller, visited, order)?;
                }
                visited[i] = 2;
                order.push(i);
                Some(())
            }
            for i in 0..num_kmers {
                dfs(i, &smaller, &mut visited, &mut order)?;
            }
            order.reverse();
        }

        let mut idx = vec![(0, Direction::Leftmost); num_kmers];
        for (i, &x) in order.iter().enumerate() {
            idx[x].0 = i;
            idx[x].1 = directions[x].unwrap_or(Direction::Leftmost);
        }

        Some(ExplicitDirectedOrder {
            k: self.k,
            sigma: self.sigma,
            idx,
        })
    }
}

/// returns number of selected positions in the text.
pub fn text_localscheme(text: &[u8], w: usize, k: usize, ls: &LocalScheme) -> usize {
    assert_eq!(w, ls.w);
    assert_eq!(k, ls.k);
    let l = k + w - 1;
    // eprintln!("k: {}", k);
    // eprintln!("w: {}", w);
    // eprintln!("l: {}", l);
    let mut taken = vec![];
    for (i, lmer) in text.windows(l).enumerate() {
        let j = ls.get(lmer);
        taken.push(i + j);
    }
    taken.sort();
    taken.dedup();
    taken.len()
}

pub fn bd_anchor(s: &[u8], r: usize) -> usize {
    assert!(r <= s.len());
    let mut best = 0;
    for i in 1..s.len().saturating_sub(r) {
        for j in 0..s.len() {
            if s[(i + j) % s.len()] != s[(best + j) % s.len()] {
                if s[(i + j) % s.len()] < s[(best + j) % s.len()] {
                    best = i;
                }
                break;
            }
        }
    }
    best
}

/// 1. Select kmers whose k0-minimizer is at their start or end.
/// 2. From these, select a minimal one.
// TODO: Same hash or distinct hashes?
pub fn miniception(s: &[u8], k: usize, k0: usize, o: &impl Order) -> usize {
    let w = s.len() - k + 1;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    let w0 = k - k0;
    let _l0 = w0 + k0 - 1;

    s.windows(k)
        .enumerate()
        .filter(|(_, kmer)| {
            let i = random_minimizer(kmer, k0);
            assert!(i <= k - k0);
            i == 0 || i == k - k0
        })
        .min_by_key(|&(_i, w)| o.key(w))
        .unwrap()
        .0
}

pub fn text_miniception<'a>(
    text: &'a [u8],
    w: usize,
    k: usize,
    k0: usize,
    o: &'a impl Order,
) -> impl Minimizers + 'a {
    // The number of k0-mers in a kmer.
    let w0 = k - k0;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    // Queue of all k0-mers.
    let mut q0 = MonotoneQueue::new();
    // Queue of filtered k-mers.
    let mut q = MonotoneQueue::new();

    // i: position of lmer
    // j: position of kmer
    // j0: position of k0mer

    // 1: init k0-mers.
    let mut k0mers = text.windows(k0).enumerate();
    for (j0, k0mer) in k0mers.by_ref().take(w0) {
        q0.push(j0, o.key(k0mer));
    }

    // 2: init k-mers.
    let mut kmers = text.windows(k).enumerate().zip(k0mers);
    for ((j, kmer), (j0, k0mer)) in kmers.by_ref().take(w - 1) {
        q0.push(j0, o.key(k0mer));
        let min_pos = q0.pop(j).unwrap().0;
        if min_pos == j || min_pos == j + w0 {
            q.push(j, o.key(kmer));
        }
    }

    // 3: Iterate l-mers.
    kmers.enumerate().map(move |(i, ((j, kmer), (j0, k0mer)))| {
        q0.push(j0, o.key(k0mer));
        let min_pos = q0.pop(j).unwrap().0;
        if min_pos == j || min_pos == j + w0 {
            q.push(j, o.key(kmer));
        }

        q.pop(i).unwrap().0
    })
}

/// Sort filtered kmers by:
/// - first prefixes by h(k0), then suffixes by h(k0)
///
/// For small k, other orderings are worse, such as:
/// - first prefixes by h(k0), then suffixes by -h(k0)
/// - first prefixes by -h(k0), then suffixes by h(k0)
/// - first prefixes by -h(k0), then suffixes by -h(k0)
///
/// Sorting by h(k0) directly is same as mod_minimizers and best for large k.
pub fn miniception_new(s: &[u8], k: usize, k0: usize, o: &impl Order) -> usize {
    let w = s.len() - k + 1;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    let w0 = k - k0;
    let _l0 = w0 + k0 - 1;

    s.windows(k)
        .enumerate()
        .filter_map(|(i, kmer)| {
            let j = random_minimizer(kmer, k0);
            assert!(j <= k - k0);
            if j == 0 || j == k - k0 {
                Some((i, kmer, j == 0, o.key(&kmer[j..j + k0])))
            } else {
                None
            }
        })
        .min_by_key(|&(_i, w, b, hk0)| (b, hk0, o.key(w)))
        .unwrap()
        .0
}

// TODO: dedup
pub fn text_miniception_new<'a>(
    text: &'a [u8],
    w: usize,
    k: usize,
    k0: usize,
    o: &'a impl Order,
) -> impl Minimizers + 'a {
    // The number of k0-mers in a kmer.
    let w0 = k - k0;
    assert!(k0 >= k.saturating_sub(w));
    assert!(k0 <= k);

    // Queue of all k0-mers.
    let mut q0 = MonotoneQueue::new();
    // Queue of filtered k-mers.
    let mut q = MonotoneQueue::new();

    // i: position of lmer
    // j: position of kmer
    // j0: position of k0mer

    // 1: init k0-mers.
    let mut k0mers = text.windows(k0).enumerate();
    for (j0, k0mer) in k0mers.by_ref().take(w0) {
        q0.push(j0, o.key(k0mer));
    }

    // 2: init k-mers.
    let mut kmers = text.windows(k).enumerate().zip(k0mers);
    for ((j, kmer), (j0, k0mer)) in kmers.by_ref().take(w - 1) {
        q0.push(j0, o.key(k0mer));
        let min_pos = q0.pop(j).unwrap().0;
        assert!(j <= min_pos && min_pos <= j + w0);
        if min_pos == j || min_pos == j + w0 {
            // TODO: Test without h(kmer)?
            q.push(
                j,
                (
                    min_pos == j,
                    o.key(&text[min_pos..min_pos + k0]),
                    o.key(kmer),
                ),
            );
        }
    }

    // 3: Iterate l-mers.
    kmers.enumerate().map(move |(i, ((j, kmer), (j0, k0mer)))| {
        q0.push(j0, o.key(k0mer));
        let min_pos = q0.pop(j).unwrap().0;
        if min_pos == j || min_pos == j + w0 {
            q.push(
                j,
                (
                    min_pos == j,
                    o.key(&text[min_pos..min_pos + k0]),
                    o.key(kmer),
                ),
            );
        }

        q.pop(i).unwrap().0
    })
}

pub fn robust_biminimizer(s: &[u8], k: usize, last: &mut usize, o: &impl Order) -> usize {
    let mut vals = s
        .windows(k)
        .enumerate()
        .map(|(i, w)| (o.key(w), i))
        .collect_vec();
    vals.sort();
    let i1 = vals[0].1;
    if vals.len() == 1 {
        *last = i1;
        return i1;
    }
    let i2 = vals[1].1;
    if *last == i1 + 1 || *last == i2 + 1 {
        *last -= 1;
    } else {
        *last = i1.max(i2);
    }
    *last
}

/// Find minimal t-mer at pos idx. Then select idx % w.
/// NOTE: THIS IS NOT A FORWARD SCHEME.
pub fn lr_minimizer(s: &[u8], k: usize, t: usize, o: &impl Order) -> usize {
    let l = s.len();
    let w = l - k + 1;
    let idx = s
        .windows(t)
        .enumerate()
        .min_by_key(|&(_i, w)| o.key(w))
        .unwrap()
        .0;
    let i = if idx >= w { idx - w } else { idx };
    assert!(i + k <= l);
    i
}

/// NOTE: THIS IS NOT A FORWARD SCHEME.
pub fn text_lr_minimizers<'a>(
    text: &'a [u8],
    w: usize,
    k: usize,
    t: usize,
    o: &'a impl Order,
) -> impl Minimizers + 'a {
    let mut q = MonotoneQueue::new();
    let l = w + k - 1;
    let wt = l - t + 1;
    let mut tmers = text.windows(t).enumerate();

    for (j, tmer) in tmers.by_ref().take(wt - 1) {
        q.push(j, o.key(tmer));
    }
    // i: position of lmer
    // j: position of tmer
    tmers.enumerate().map(move |(i, (j, tmer))| {
        q.push(j, o.key(tmer));
        let idx = q.pop(i).unwrap().0 - i;
        i + if idx >= w { idx - w } else { idx }
    })
}

/// Find minimal t-mer at pos idx. Then select idx % w.
/// NOTE: THIS IS NOT A FORWARD SCHEME.
pub fn mod_minimizer(s: &[u8], k: usize, t: usize, o: &impl Order) -> usize {
    let l = s.len();
    let w = l - k + 1;
    let idx = s
        .windows(t)
        .enumerate()
        .min_by_key(|&(_i, w)| o.key(w))
        .unwrap()
        .0;
    idx % w
}

/// FastMod32, using the low 32 bits of the hash.
/// Taken from https://github.com/lemire/fastmod/blob/master/include/fastmod.h
#[derive(Copy, Clone, Debug)]
pub struct FM32 {
    d: u64,
    m: u64,
}
impl FM32 {
    fn new(d: usize) -> Self {
        assert!(d <= u32::MAX as usize);
        Self {
            d: d as u64,
            m: (u64::MAX / d as u64).wrapping_add(1),
        }
    }
    fn reduce(self, h: usize) -> usize {
        let lowbits = self.m.wrapping_mul(h as u64);
        ((lowbits as u128 * self.d as u128) >> 64) as usize
    }
}

/// NOTE: THIS IS NOT A FORWARD SCHEME.
pub fn text_mod_minimizers<'a>(
    text: &'a [u8],
    w: usize,
    k: usize,
    t: usize,
    o: &'a impl Order,
) -> impl Minimizers + 'a {
    let mut q = MonotoneQueue::new();
    let l = w + k - 1;
    let wt = l - t + 1;
    let mut tmers = text.windows(t).enumerate();

    let fastmod_w = FM32::new(w);

    for (j, tmer) in tmers.by_ref().take(wt - 1) {
        q.push(j, o.key(tmer));
    }
    // i: position of lmer
    // j: position of tmer
    tmers.enumerate().map(move |(i, (j, tmer))| {
        q.push(j, o.key(tmer));
        i + fastmod_w.reduce(q.pop(i).unwrap().0 - i)
    })
}

/// Asymptotic (in k) optimal minimizers:
/// - Assume k=x*w; sum all i mod w positions.
/// - Take minimizer if first coordinate is max or at most sigma away from max.
/// - If multiple, take random.
///
/// - Ignore for k < w.
/// - Ignore trailing chars.
///
/// Simplified:
/// - For all k-mers compute the rotational sum of k/w w-strided values.
/// - Take k-mer that maximizes the first coordinate.
pub fn rot_minimizer(s: &[u8], k: usize) -> usize {
    let w = s.len() - k + 1;
    assert!(k >= w);
    (0..=s.len() - k)
        .max_by_key(|&i| {
            s[i..i + k]
                .iter()
                .step_by(w)
                .map(|&x| x as u64)
                .sum::<u64>()
        })
        .unwrap()
}

pub fn decycling_minimizer_init(k: usize) -> Vec<f32> {
    // (0..k).map(|i| (i as f32).sin()).collect()
    (0..k)
        .map(|i| (i as f32 / k as f32 * 2. * PI).sin())
        .collect()
}

/// Check if this is in the decyling set.
/// Note that we ignore the edge case where ix=0 and ix2=0.
pub fn is_in_decycling_set(kmer: &[u8], cs: &Vec<f32>) -> bool {
    let ix = zip(kmer, cs).map(|(&x, c)| x as f32 * c).sum::<f32>();
    let ix2 = *kmer.last().unwrap() as f32 * cs[0]
        + zip(kmer, &cs[1..]).map(|(&x, c)| x as f32 * c).sum::<f32>();
    ix > 0. && ix2 <= 0.
}
pub fn is_in_double_decycling_set(kmer: &[u8], cs: &Vec<f32>) -> bool {
    let ix = zip(kmer, cs).map(|(&x, c)| x as f32 * c).sum::<f32>();
    let ix2 = *kmer.last().unwrap() as f32 * cs[0]
        + zip(kmer, &cs[1..]).map(|(&x, c)| x as f32 * c).sum::<f32>();
    ix < 0. && ix2 >= 0.
}

/// Decycling minimizer using Mykkelvelt embedding
pub fn decycling_minimizer(s: &[u8], k: usize, cs: &Vec<f32>, o: &impl Order) -> usize {
    s.windows(k)
        .enumerate()
        .min_by_key(|&(_i, kmer)| (!is_in_decycling_set(kmer, cs), o.key(kmer)))
        .unwrap()
        .0
}

/// Decycling minimizer using Mykkelvelt embedding
pub fn double_decycling_minimizer(s: &[u8], k: usize, cs: &Vec<f32>, o: &impl Order) -> usize {
    s.windows(k)
        .enumerate()
        .min_by_key(|&(_i, kmer)| {
            (
                !is_in_decycling_set(kmer, cs),
                !is_in_double_decycling_set(kmer, cs),
                o.key(kmer),
            )
        })
        .unwrap()
        .0
}

/// Finds the best order for a minimizer scheme.
///
/// Returns (# selected, # total), order.
pub fn bruteforce_minimizer(k: usize, w: usize, sigma: usize) -> ((usize, usize), ExplicitOrder) {
    let (text, len) = cyclic_exact_density_string(k, w, sigma, true);
    let num_kmers = sigma.pow(k as u32);

    let mut perms = 1;
    for i in 1..=num_kmers {
        perms *= i;
    }
    eprintln!("Num permutations: {num_kmers}! = {perms}");

    let best = (0..num_kmers)
        .permutations(num_kmers)
        .par_bridge()
        .map(|perm| {
            let o = ExplicitOrder {
                k,
                sigma,
                idx: perm,
            };
            let cnt = cyclic_minimizers(&text, len, w, k, &o);
            (cnt, o)
        })
        .min_by_key(|x| x.0)
        .unwrap();
    let (cnt, o) = best;
    ((cnt, len), o)
}

/// Finds the best directed minimizer scheme.
///
/// Returns (# selected, # total), order.
pub fn bruteforce_directed_minimizer(
    k: usize,
    w: usize,
    sigma: usize,
) -> ((usize, usize), ExplicitDirectedOrder) {
    let (text, len) = cyclic_exact_density_string(k, w, sigma, true);
    let num_kmers = sigma.pow(k as u32);

    let mut perms = 1;
    for i in 1..=num_kmers {
        perms *= i * 2;
    }
    eprintln!("Num permutations: {num_kmers}! * 2^{num_kmers} = {perms}");

    let best = (0..num_kmers)
        .permutations(num_kmers)
        .par_bridge()
        .map(|perm| {
            repeat_n([Direction::Leftmost, Direction::Rightmost], num_kmers)
                .multi_cartesian_product()
                .map(|directions| {
                    let o = ExplicitDirectedOrder {
                        k,
                        sigma,
                        idx: zip(perm.iter().copied(), directions.iter().copied()).collect_vec(),
                    };
                    let cnt = cyclic_minimizers(&text, len, w, k, &o);
                    (cnt, o)
                })
                .min_by_key(|x| x.0)
                .unwrap()
        })
        .min_by_key(|x| x.0)
        .unwrap();
    let (cnt, o) = best;
    ((cnt, len), o)
}

/// Finds the best local scheme by trying all mappings sigma^l -> [w].
///
/// Returns (# selected, # total), mapping.
pub fn bruteforce_local_scheme(k: usize, w: usize, sigma: usize) -> ((usize, usize), LocalScheme) {
    let l = k + w - 1;
    let num_lmers = sigma.pow(l as u32);
    let text = exact_density_string(k, w, sigma, false);
    eprintln!("Num lmers: {}", num_lmers);
    eprintln!("Num maps : {}", w.pow(num_lmers as u32));
    eprintln!("text len : {}", text.len());
    let best = repeat_n(0..w as u8, num_lmers)
        .multi_cartesian_product()
        .par_bridge()
        .map(|map| {
            let ls = LocalScheme { k, w, sigma, map };
            let cnt = text_localscheme(&text, w, k, &ls);
            (cnt, ls)
        })
        .min_by_key(|x| x.0)
        .unwrap();
    let (cnt, m) = best;
    ((cnt, text.windows(k).len()), m)
}

pub struct SuperKmer {
    /// Start position in text of the first basepair.
    pub start: usize,
    /// Past-the-end index of the start of the last kmer.
    pub kmer_end: usize,
    /// Past-the-end index of the full super-kmer.
    pub super_kmer_end: usize,
    /// The position of the minimizer in the super-kmer.
    pub minimizer_pos: usize,
}

/// Takes an iterator over the minimize positions in a string and returns an
/// iterator over super kmers, i.e. substrings sharing the same minimizer.
///
/// Returns (start, end, &[u8])
pub fn super_kmers<'a>(
    text: &'a [u8],
    k: usize,
    minimizer_positions: impl Minimizers + 'a,
) -> impl Iterator<Item = SuperKmer> + 'a {
    let mut first = None;
    minimizer_positions.enumerate().filter_map(move |(i, pos)| {
        let Some((first_i, first_pos)) = first else {
            // Initial set.
            first = Some((i, pos));
            return None;
        };
        if i + 1 == text.len() {
            // Save last.
            return Some(SuperKmer {
                start: first_i,
                kmer_end: i + 1,
                super_kmer_end: i + 1 + k,
                minimizer_pos: first_pos,
            });
        }
        // Extend?
        if pos == first_pos {
            return None;
        } else {
            // Start new.
            first = Some((i, pos));
            // Save previous.
            return Some(SuperKmer {
                start: first_i,
                kmer_end: i,
                super_kmer_end: i + k,
                minimizer_pos: first_pos,
            });
        }
    })
}
