#![feature(exclusive_range_pattern)]
use std::io::{stdout, Write};

use clap::Parser;
use colored::Colorize;
use itertools::Itertools;

/// Generate a random string.
fn generate_random_string(n: usize, sigma: usize) -> Vec<u8> {
    (0..n)
        .map(|_| (rand::random::<usize>() % sigma) as u8)
        .collect()
}

fn minimal_rotation_pos(s: &[u8], r: usize) -> usize {
    let mut best = 0;
    for i in 1..s.len() - r {
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

fn minimizer_pos(s: &[u8], k: usize) -> usize {
    s.windows(k)
        .enumerate()
        .min_by_key(|(_, w)| fxhash::hash64(w))
        .unwrap()
        .0
}

fn biminimizer_pos(s: &[u8], k: usize, last: &mut usize) -> usize {
    let h1 = |x: &[u8]| fxhash::hash64(x);
    let h2 = |x: &[u8]| fxhash::hash64(&fxhash::hash64(x));
    let i1 = s
        .windows(k)
        .enumerate()
        .min_by_key(|(_, w)| h1(w))
        .unwrap()
        .0;
    let i2 = s
        .windows(k)
        .enumerate()
        .min_by_key(|(_, w)| h2(w))
        .unwrap()
        .0;
    if *last == i1 + 1 || *last == i2 + 1 {
        *last -= 1;
    } else {
        *last = i1.max(i2);
    }
    *last
}

fn reduced_minimizer_pos(s: &[u8], k: usize, r: usize) -> usize {
    let idx = s
        .windows(r)
        .enumerate()
        .min_by_key(|(_, w)| fxhash::hash64(w))
        .unwrap()
        .0;
    let a = (2 * s.len() + r) * 1 / 3;
    let b = (2 * s.len() + r) * 2 / 3;
    let pos = 2 * idx + r;
    let i = if pos < a {
        idx
    } else if pos < b {
        idx - (k - r) / 2
    } else {
        idx - (k - r)
    };
    // let i = if idx > s.len() - k {

    //     // if idx >= k - r {
    //     //     idx - k + r
    //     // } else {
    //     //     idx - k / 2 + r
    //     // }
    // } else {
    //     idx
    // };
    assert!(i <= s.len());
    assert!(
        i + k <= s.len(),
        "len={} idx={idx} k={k} r={r} pos={pos} i={i}",
        s.len()
    );
    i
}

#[derive(clap::ValueEnum, Clone, Copy, Default)]
enum Type {
    Minimizer,
    ReducedMinimizer,
    BiMinimizer,
    #[default]
    BdAnchor,
}

/// Print expected and measured densities for different l and r.
#[derive(clap::Parser)]
struct Args {
    /// Length of the generated random string.
    #[arg(short, default_value_t = 100000)]
    n: usize,
    /// Alphabet size.
    #[arg(short, long, default_value_t = 4)]
    sigma: usize,
    /// Test kmers up to l.
    #[arg(short, default_value_t = 40)]
    l: usize,
    /// Test r (reduction constant) up to this.
    #[arg(short, default_value_t = 10)]
    r: usize,
    /// Minimizer length.
    #[arg(short, default_value_t = 10)]
    k: usize,
    /// The type of anchor to use.
    #[arg(short, value_enum, default_value_t = Type::default())]
    t: Type,
}

fn main() {
    let args = Args::parse();
    let t = generate_random_string(args.n, args.sigma);

    match args.t {
        Type::BdAnchor => {
            print!(" l ");
            for r in args.r..=args.r {
                print!("        r={r:>2}");
            }
            println!();
            for l in args.l..=args.l {
                print!("{l:>2}:");
                let mut last_density = 2.;
                for r in args.r..=args.r.min(l - 1) {
                    let anchors: Vec<_> = t
                        .windows(l)
                        .enumerate()
                        .map(|(i, w)| i + minimal_rotation_pos(w, r))
                        .dedup()
                        .collect();
                    let exp = 2. / (l + 1 - r) as f64;
                    let density = anchors.len() as f64 / (t.len() - l + 1) as f64;
                    let density_str = if density < last_density {
                        format!("{density:>.3}").green()
                    } else {
                        format!("{density:>.3}").red()
                    };
                    print!(" {}/{density_str}", format!("{exp:>.3}").dimmed(),);
                    last_density = density;
                    stdout().flush().unwrap();
                }
                println!();
            }
        }

        Type::Minimizer => {
            print!(" l ");
            for k in args.k..=args.k {
                print!("        k={k:>2}");
            }
            println!();
            for l in args.l..=args.l {
                print!("{l:>2}:");
                let mut last_density = 2.;
                for k in args.k..=args.k {
                    let anchors: Vec<_> = t
                        .windows(l)
                        .enumerate()
                        .map(|(i, w)| i + minimizer_pos(w, k))
                        .dedup()
                        .collect();
                    let exp = 2. / (l + 2 - k) as f64;
                    let density = anchors.len() as f64 / (t.len() - l + 1) as f64;
                    let density_str = if density < last_density {
                        format!("{density:>.3}").green()
                    } else {
                        format!("{density:>.3}").red()
                    };
                    print!(" {}/{density_str}", format!("{exp:>.3}").dimmed(),);
                    last_density = density;
                    stdout().flush().unwrap();
                }
                println!();
            }
        }

        Type::BiMinimizer => {
            print!(" l ");
            for k in args.k..=args.k {
                print!("        k={k:>2}");
            }
            println!();
            for l in args.l..=args.l {
                print!("{l:>2}:");
                let mut last_density = 2.;
                for k in args.k..=args.k {
                    let last = &mut 0;
                    let anchors: Vec<_> = t
                        .windows(l)
                        .enumerate()
                        .map(|(i, w)| i + biminimizer_pos(w, k, last))
                        .dedup()
                        .collect();
                    let exp = 2. / (l + 2 - k) as f64;
                    let density = anchors.len() as f64 / (t.len() - l + 1) as f64;
                    let density_str = if density < last_density {
                        format!("{density:>.3}").green()
                    } else {
                        format!("{density:>.3}").red()
                    };
                    print!(" {}/{density_str}", format!("{exp:>.3}").dimmed(),);
                    last_density = density;
                    stdout().flush().unwrap();
                }
                println!();
            }
        }

        Type::ReducedMinimizer => {
            print!(" l ");
            for k in args.k..=args.k {
                print!("        k={k:>2}");
            }
            println!();
            for l in args.l..=args.l {
                print!("{l:>2}:");
                let mut last_density = 2.;
                for k in args.k..=args.k {
                    let anchors: Vec<_> = t
                        .windows(l)
                        .enumerate()
                        .map(|(i, w)| i + reduced_minimizer_pos(w, k, args.r))
                        .dedup()
                        .collect();
                    let exp = 2. / (l + 2 - k) as f64;
                    let density = anchors.len() as f64 / (t.len() - l + 1) as f64;
                    let density_str = if density < last_density {
                        format!("{density:>.3}").green()
                    } else {
                        format!("{density:>.3}").red()
                    };
                    print!(" {}/{density_str}", format!("{exp:>.3}").dimmed(),);
                    last_density = density;
                    stdout().flush().unwrap();
                }
                println!();
            }
        }
    }
}
