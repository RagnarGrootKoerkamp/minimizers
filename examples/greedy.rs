//! Manual implementations of some greedymini schemes.

use std::cmp::Reverse;

use clap::Parser;
use itertools::Itertools;
use minimizers::{
    collect_stats, generate_random_string,
    schemes::{greedy::GreedyP, *},
    Params,
};

#[derive(clap::Parser)]
struct Args {
    #[clap(short)]
    k: usize,
    #[clap(short)]
    w: usize,
}

fn main() {
    let Args { k, w } = Args::parse();
    let sigma = 4;

    let m = M(GreedyP).build(w, k, sigma);
    let text = &generate_random_string(10000000, sigma);
    let (d, _pos, _dist, _trans, counts, dist_counts) = collect_stats(w, k, text, &*m);
    let lb = (w + k).div_ceil(w) as f64 / (w + k) as f64;
    let k2 = w.next_multiple_of(k) + 1;
    let lb2 = (w + k2).div_ceil(w) as f64 / (w + k2) as f64;
    let lb = lb.max(lb2);
    eprintln!("w: {w} k: {k}");
    eprintln!("density  {d:.3}   lb: {lb:.3}",);
    eprintln!("avg dist {:.3}   ub: {:.3}", 1. / d, 1. / lb,);
    // eprintln!("pos: {:?}", pos);
    // eprintln!("dist: {:?}", &dist[w + 1..]);
    eprintln!("kmer   count  |     1   2   3  dist(%)");
    let mut items = counts.iter().collect_vec();
    items.sort_by_key(|x| Reverse(x.1));
    for (kmer, c) in items {
        eprint!("{kmer:>0k$b}: {:>7}  |    ", c / 1000);
        let dist_counts = dist_counts.get(kmer).unwrap();
        for i in 1..=w {
            let c = *dist_counts.get(i).unwrap() as f32 / *c as f32 * 100.0;
            eprint!("{c:>2.0}  ");
        }
        eprintln!();
    }
}
