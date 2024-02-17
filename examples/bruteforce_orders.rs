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
use minimizers::bruteforce_minimizer;

#[derive(clap::Parser)]
struct Args {
    #[clap(short, default_value_t = 2)]
    k: usize,
    #[clap(short, default_value_t = 2)]
    w: usize,
    #[clap(short, default_value_t = 2)]
    sigma: usize,
}

fn main() {
    let Args { k, w, sigma } = Args::parse();

    for w in 2..=w {
        eprintln!("w = {w}");
        let best = bruteforce_minimizer(k, w, sigma);

        eprintln!("  best cnt : {} / {}", best.0 .0, best.0 .1);
        eprintln!("  density  : {}", (best.0 .0) as f32 / (best.0 .1) as f32);
        eprintln!("  best perm: {:?}", best.1.idx);
    }
}
