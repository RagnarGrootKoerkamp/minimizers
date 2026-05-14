//! Manual implementations of some greedymini schemes.

use std::cmp::Reverse;

use clap::Parser;
use colored::Colorize;
use itertools::Itertools;
use minimizers::{
    collect_stats, generate_random_string,
    order::Order,
    schemes::{greedy::GreedyP, *},
    Params, SamplingScheme,
};

#[derive(clap::Parser)]
struct Args {
    #[clap(short)]
    k: usize,
    #[clap(short)]
    w: usize,

    /// Input size x * 10MB
    #[clap(short, default_value_t = 10)]
    x: usize,

    #[clap(long)]
    cycle_stats: Option<usize>,

    #[clap(long)]
    tmers: Vec<String>,
}

fn main() {
    let Args {
        k,
        w,
        x,
        tmers,
        cycle_stats,
    } = Args::parse();
    let n = x * 1000000;
    let sigma = 2;

    let text = &generate_random_string(n, sigma);

    let m: Box<dyn SamplingScheme>;
    if !tmers.is_empty() {
        eprintln!("tmers: {:?}", tmers);
        let order = F {
            tmers: tmers
                .iter()
                .map(|x| {
                    let (dir, idx) = match x.as_bytes()[0] {
                        b'_' => (Dir::StartEnd, 1),
                        b'<' => (Dir::Start, 1),
                        b'>' => (Dir::End, 1),
                        _ => (Dir::EndStart, 0),
                    };
                    return (
                        dir,
                        x[idx..].as_bytes().iter().map(|x| x - b'0').collect_vec(),
                    );
                })
                .collect(),
        };
        // eprintln!("order: {:?}", order);
        m = Box::new(Minimizer::new(k, w, order));
    } else {
        m = M(GreedyP).build(w, k, sigma);
    }
    let (d, _pos, _dist, _trans, counts, dist_counts) = collect_stats(w, k, text, &*m);
    let lb = (w + k).div_ceil(w) as f64 / (w + k) as f64;
    let k2 = k.next_multiple_of(w) + 1;
    let lb2 = (w + k2).div_ceil(w) as f64 / (w + k2) as f64;
    let lb = lb.max(lb2);
    eprintln!("w: {w} k: {k}");
    eprintln!("density  {d:.3}   lb: {lb:.3}",);
    eprintln!("avg dist {:.3}   ub: {:.3}", 1. / d, 1. / lb,);
    eprintln!();
    eprint!(
        "{:<k$}        {:<y$} dist_from (%) {:<x$} dist_to (%){:<y$}avg-from   avg-to     avg    {:<y$} pos (%)",
        "",
        "",
        "",
        "",
        "",
        y = 2 * w - 3,
        x = 4 * w - 7
    );
    eprintln!();
    eprint!("{:<k$}  count  |  ", "kmer", k = k + 1);
    for i in 1..=w {
        eprint!("{i:>3} ");
    }
    eprint!("  |  ");
    for i in 1..=w {
        eprint!("{i:>3} ");
    }
    eprint!("  |  ");
    eprint!("                         ");
    eprint!("  |  ");
    for i in 0..w {
        eprint!("{i:>3} ");
    }
    eprintln!();
    let mut items = counts.iter().map(|(&k, &c)| (k, c)).collect_vec();

    fn avg(d: &Vec<usize>) -> f32 {
        let sum = d.iter().sum::<usize>();
        let mut avg = 0;
        for (i, &v) in d.iter().enumerate() {
            avg += i * v;
        }
        avg as f32 / sum as f32
    }

    items.sort_by_key(|x| Reverse(x.1));
    items.sort_by(|(kmera, _cnta), (kmerb, _cntb)| {
        let (d_to, d_from, _) = dist_counts.get(&kmera).unwrap();
        let a1 = (avg(d_to) + avg(d_from)) / 2.0;
        let (d_to, d_from, _) = dist_counts.get(&kmerb).unwrap();
        let a2 = (avg(d_to) + avg(d_from)) / 2.0;
        a2.total_cmp(&a1)
    });

    fn format_dist(c: usize, v: usize, m: usize) {
        let c = v as f32 / c as f32 * 100.0;
        if c < 1.0 {
            eprint!("{} ", format!("{c:>3.0}").truecolor(20, 20, 20));
        } else if c < 10.0 {
            eprint!("{} ", format!("{c:>3.0}").truecolor(130, 130, 130));
        } else if v == m {
            eprint!("{} ", format!("{c:>3.0}").bold().green());
        } else {
            eprint!("{c:>3.0} ");
        }
    }

    for (kmer, c) in items {
        eprint!("{kmer:>0k$b} {:>7}  |  ", c / 1000);
        let (dist_to_counts, dist_from_counts, pos_counts) = dist_counts.get(&kmer).unwrap();
        let m = *dist_to_counts.iter().max().unwrap();
        let mut avg_dist_to = 0;
        for (i, &v) in dist_to_counts.iter().enumerate().skip(1) {
            format_dist(c, v, m);
            avg_dist_to += i * v;
        }
        eprint!("  |  ");
        let m = *dist_from_counts.iter().max().unwrap();
        let mut avg_dist_from = 0;
        for (i, &v) in dist_from_counts.iter().enumerate().skip(1) {
            format_dist(c, v, m);
            avg_dist_from += i * v;
        }
        eprint!("  |  ");
        eprint!("{:>5.3}", avg_dist_to as f32 / c as f32);
        eprint!("  |  ");
        eprint!("{:>5.3}", avg_dist_from as f32 / c as f32);
        eprint!("  |  ");
        eprint!(
            "{:>5.3}",
            (avg_dist_from + avg_dist_to) as f32 / (2 * c) as f32
        );
        eprint!("  |  ");
        let sum_pos = pos_counts.iter().sum::<usize>();
        for i in 0..w {
            let c = pos_counts[i] as f32 / sum_pos as f32 * 100.0;
            eprint!("{c:>3.0} ");
        }
        eprintln!();
    }

    if let Some(l) = cycle_stats {
        let (density, stats) = minimizers::cycle_stats(l, w, text, &*m);
        eprintln!("Cycle density: {density}");
        eprintln!("Samples:");
        for (i, c) in stats.iter().enumerate() {
            eprintln!("{i:>2}: {c:>6}");
        }
    }
}

#[derive(Debug)]
enum Dir {
    StartEnd,
    EndStart,
    Start,
    End,
}

#[derive(Debug)]
struct F {
    tmers: Vec<(Dir, Vec<u8>)>,
}

impl Order for F {
    /// The matching tmer, and dist to next/prev occ of it.
    type T = (usize, usize, usize);

    fn key(&self, kmer: &[u8]) -> Self::T {
        let find_start = |tmer| -> Option<usize> {
            kmer.starts_with(tmer).then(|| {
                // next occ of tmer
                kmer[1..]
                    .windows(tmer.len())
                    .position(|w| w == tmer)
                    .map_or(1000, |i| 100 - i)
            })
        };
        let find_end = |tmer| -> Option<usize> {
            kmer.ends_with(tmer).then(|| {
                // next occ of tmer
                kmer[..kmer.len() - 2]
                    .windows(tmer.len())
                    .rposition(|w| w == tmer)
                    .map_or(1000, |i| i + 1)
            })
        };

        for (i, (dir, tmer)) in self.tmers.iter().enumerate() {
            match *dir {
                Dir::StartEnd => {
                    if let Some(next) = find_start(tmer) {
                        return (i, 0, next);
                    }
                    if let Some(next) = find_end(tmer) {
                        return (i, 1, next);
                    }
                }
                Dir::EndStart => {
                    if let Some(next) = find_end(tmer) {
                        return (i, 0, next);
                    }
                    if let Some(next) = find_start(tmer) {
                        return (i, 1, next);
                    }
                }
                Dir::Start => {
                    if let Some(next) = find_start(tmer) {
                        return (i, 0, next);
                    }
                }
                Dir::End => {
                    if let Some(next) = find_end(tmer) {
                        return (i, 0, next);
                    }
                }
            }
        }
        return (1000, 0, 0);
    }
}
