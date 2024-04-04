//! Some code to bruteforce over all possible orders for minimizers.
//!
//! Results:
//! w   k   s   density
//! 2   2   2   12/18   = 0.666
//! 2   2   3   53/83   = 0.638
//! 2   3   2   21/34   = 0.617
//! 3   2   2   17/35   = 0.485
//! 3   2   3   114/246 = 0.463
//! 3   3   2   30/67   = 0.447
//! 4   2   2   
//! 4   2   3   114/246 = 0.463
//! 4   3   2   30/67   = 0.447

use clap::Parser;
use minimizers::{
    bruteforce_directed_minimizer, bruteforce_local_scheme, bruteforce_minimizer, ilp_scheme,
    LocalScheme,
};

#[derive(clap::ValueEnum, Clone)]
enum Method {
    BFMinimizer,
    BFDirectedMinimizer,
    BFLocalScheme,
    ILP,
    Verify,
}

#[derive(clap::Parser)]
struct Args {
    #[clap(value_enum)]
    method: Method,
    #[clap(short, default_value_t = 2)]
    k: usize,
    #[clap(short, default_value_t = 2)]
    w: usize,
    #[clap(short, default_value_t = 2)]
    sigma: usize,
    #[clap(long)]
    forward: bool,
}

fn main() {
    let Args {
        k,
        w,
        sigma,
        method,
        forward,
    } = Args::parse();

    match method {
        Method::BFMinimizer => {
            let best = bruteforce_minimizer(k, w, sigma);
            eprintln!("  best cnt : {} / {}", best.0 .0, best.0 .1);
            eprintln!("  density  : {}", (best.0 .0) as f32 / (best.0 .1) as f32);
            eprintln!("  best perm: {:?}", best.1.idx);
        }
        Method::BFDirectedMinimizer => {
            let best = bruteforce_directed_minimizer(k, w, sigma);
            eprintln!("  best cnt : {} / {}", best.0 .0, best.0 .1);
            eprintln!("  density  : {}", (best.0 .0) as f32 / (best.0 .1) as f32);
            eprintln!("  best perm: {:?}", best.1.idx);
        }
        Method::BFLocalScheme => {
            let best = bruteforce_local_scheme(k, w, sigma);
            eprintln!("  best cnt : {} / {}", best.0 .0, best.0 .1);
            eprintln!("  density  : {}", (best.0 .0) as f32 / (best.0 .1) as f32);
            eprintln!("  best map: {:?}", best.1.map);
        }
        Method::ILP => {
            let best = ilp_scheme::best_local_scheme(k, w, sigma, forward);
            eprintln!("  best cnt : {} / {}", best.0 .0, best.0 .1);
            eprintln!("  density  : {}", (best.0 .0) as f32 / (best.0 .1) as f32);
            eprintln!("  best map: {:?}", best.1.map);
        }
        Method::Verify => {
            let best = vec![
                LocalScheme {
                    k: 2,
                    w: 2,
                    sigma: 2,
                    map: vec![0, 1, 0, 0, 1, 1, 0, 1],
                },
                LocalScheme {
                    k: 2,
                    w: 2,
                    sigma: 3,
                    map: vec![
                        1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0,
                        0, 1,
                    ],
                },
                LocalScheme {
                    k: 2,
                    w: 2,
                    sigma: 4,
                    map: vec![
                        1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1,
                        1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1,
                        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1,
                    ],
                },
                LocalScheme {
                    k: 2,
                    w: 3,
                    sigma: 2,
                    map: vec![0, 2, 1, 1, 0, 2, 0, 0, 1, 2, 1, 1, 2, 2, 1, 2],
                },
                LocalScheme {
                    k: 2,
                    w: 3,
                    sigma: 2,
                    map: vec![0, 2, 1, 1, 0, 0, 0, 0, 1, 2, 1, 1, 2, 2, 0, 1],
                },
                LocalScheme {
                    k: 2,
                    w: 4,
                    sigma: 2,
                    map: vec![
                        0, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 2, 2, 1, 1, 1, 1, 2,
                        3, 2, 2, 0, 3, 1, 2,
                    ],
                },
                LocalScheme {
                    k: 2,
                    w: 4,
                    sigma: 2,
                    map: vec![
                        3, 2, 3, 1, 2, 2, 3, 3, 1, 1, 3, 1, 2, 2, 3, 2, 0, 0, 0, 0, 0, 2, 3, 0, 1,
                        1, 3, 1, 2, 2, 3, 1,
                    ],
                },
                LocalScheme {
                    k: 2,
                    w: 3,
                    sigma: 3,
                    map: vec![
                        2, 1, 2, 0, 0, 2, 1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 2, 0, 2, 1, 1, 1, 0, 0, 2, 0, 0,
                        2, 0, 0, 0, 2, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 2, 0, 2,
                        2, 0, 2, 1, 1, 2,
                    ],
                },
                LocalScheme {
                    k: 2,
                    w: 3,
                    sigma: 3,
                    map: vec![
                        2, 1, 2, 2, 0, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 0, 0, 2, 2, 0, 2, 0,
                        0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 1, 1, 1, 2, 0, 2, 1, 1, 1, 0, 0, 2, 0, 0,
                        0, 0, 0, 0, 2, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 2, 0, 2,
                        2, 2, 2, 1, 1, 2,
                    ],
                },
            ];

            for ls in best {
                eprintln!("Verifying {:?}", ls);
                if let Some(_) = ls.to_directed_order() {
                    eprintln!("  OK");
                } else {
                    eprintln!("  Failed");
                }
            }
        }
    }
}
