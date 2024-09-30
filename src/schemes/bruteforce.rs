use super::ExplicitLocalScheme;
use super::*;
use crate::de_bruijn_seq::cyclic_exact_density_string;
use crate::de_bruijn_seq::exact_density_string;
use crate::Direction;
use crate::ExplicitDirectedOrder;
use crate::ExplicitOrder;
use crate::SamplingScheme;

use itertools::repeat_n;
use itertools::Itertools;
use minimizer::Minimizer;
use rayon::iter::{ParallelBridge, ParallelIterator};
use std::iter::zip;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BruteforceP;

#[typetag::serde]
impl Params for BruteforceP {
    fn build(&self, w: usize, k: usize, sigma: usize) -> Box<dyn SamplingScheme> {
        let m = bruteforce::bruteforce_minimizer(k, w, sigma).1;
        Box::new(m)
    }
}

/// Finds the best order for a minimizer scheme.
///
/// Returns (# selected, # total), order.
pub fn bruteforce_minimizer(
    k: usize,
    w: usize,
    sigma: usize,
) -> ((usize, usize), Minimizer<ExplicitOrder>) {
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
            let minimizer = Minimizer::new(k, w, o);
            let cnt = minimizer.cyclic_text_density(&text, len);
            (cnt, minimizer)
        })
        .min_by_key(|x| x.0)
        .unwrap();
    let (cnt, minimizer) = best;
    ((cnt, len), minimizer)
}

/// Finds the best directed minimizer scheme.
///
/// Returns (# selected, # total), order.
pub fn bruteforce_directed_minimizer(
    k: usize,
    w: usize,
    sigma: usize,
) -> ((usize, usize), Minimizer<ExplicitDirectedOrder>) {
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
                    let minimizer = Minimizer::new(k, w, o);
                    let cnt = minimizer.cyclic_text_density(&text, len);
                    (cnt, minimizer)
                })
                .min_by_key(|x| x.0)
                .unwrap()
        })
        .min_by_key(|x| x.0)
        .unwrap();
    let (cnt, minimizer) = best;
    ((cnt, len), minimizer)
}

/// Finds the best local scheme by trying all mappings sigma^l -> [w].
///
/// Returns (# selected, # total), mapping.
pub fn bruteforce_local_scheme(
    k: usize,
    w: usize,
    sigma: usize,
) -> ((usize, usize), ExplicitLocalScheme) {
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
            let ls = ExplicitLocalScheme { k, w, sigma, map };
            let cnt = ls.sample(&text);
            (cnt, ls)
        })
        .min_by_key(|x| x.0)
        .unwrap();
    let (cnt, m) = best;
    ((cnt, text.windows(k).len()), m)
}
