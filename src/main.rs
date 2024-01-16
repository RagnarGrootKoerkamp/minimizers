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

/// Print expected and measured densities for different l and r.
#[derive(clap::Parser)]
struct Args {
    /// Length of the generated random string.
    #[arg(short, long, default_value = "100000")]
    n: usize,
    /// Alphabet size.
    #[arg(short, long, default_value = "4")]
    sigma: usize,
    /// Test kmers up to l.
    #[arg(short, long, default_value = "40")]
    l: usize,
    /// Test r (reduction constant) up to this.
    #[arg(short, long, default_value = "10")]
    r: usize,
}

fn main() {
    let args = Args::parse();
    let t = generate_random_string(args.n, args.sigma);
    print!(" l ");
    for r in 0..=args.r {
        print!("        r={r:>2}");
    }
    println!();
    'l: for l in 1..=args.l {
        print!("{l:>2}:");
        let mut last_density = 2.;
        for r in 0..=args.r {
            if r >= l {
                println!();
                continue 'l;
            }
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
