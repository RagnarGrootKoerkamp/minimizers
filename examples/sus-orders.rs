use clap::Parser;
use minimizers::{
    de_bruijn_seq::cyclic_exact_density_string, generate_true_random_string, order::RandomLexO,
    schemes::SusAnchorS, SamplingScheme,
};
use rand::random;

#[derive(clap::Parser)]
struct Args {
    w: usize,
    #[clap(long)]
    exact: bool,
    #[clap(long)]
    seed: Vec<usize>,
    #[clap(long)]
    v2: bool,
}

fn main() {
    let k = 1usize;
    // read from command line args
    let args = Args::parse();
    let w = args.w;
    let exact = args.exact;
    let sigma = 2;

    let len;
    let t1;
    let t2;
    if exact {
        (t1, len) = cyclic_exact_density_string(k, w, sigma, true);
        t2 = vec![];
    } else {
        len = 100000usize;
        t1 = generate_true_random_string(len, sigma);
        t2 = generate_true_random_string(len, sigma);
    }

    if !args.seed.is_empty() {
        for w in 2..=w {
            eprint!("w {w:>3}:  ");
            for &seed in &args.seed {
                let s = SusAnchorS::new(w, k, RandomLexO(seed as u128));
                let lb = (len * (w + k).div_ceil(w) / (w + k)) as isize;
                let cnt;
                if exact {
                    cnt = s.cyclic_text_samples(&t1, len).len() as isize;
                } else {
                    let mut pos = s.stream(&t1);
                    pos.dedup();
                    cnt = pos.len() as isize;
                }
                let lbd = lb as f32 / len as f32;
                let cntd = cnt as f32 / len as f32;
                eprint!("  {:>6}  {:>.6}", cnt - lb, cntd - lbd);
            }
            eprintln!();
        }
        return;
    }

    let lb = (len * (w + k).div_ceil(w) / (w + k)) as isize;

    let mut best = (isize::MAX, 0, 0);

    for i in 0..1000 {
        // let s = s.build(w, k, sigma);
        let seed = random::<u32>() as usize;
        let s = Box::new(SusAnchorS::new(w, k, RandomLexO(seed as u128)));
        let cnt;
        let cnt2;
        if exact {
            cnt = s.cyclic_text_samples(&t1, len).len() as isize;
            cnt2 = cnt;
        } else {
            let mut pos = s.stream(&t1);
            pos.dedup();
            cnt = pos.len() as isize;
            let mut pos = s.stream(&t2);
            pos.dedup();
            cnt2 = pos.len() as isize;
        }
        if cnt <= best.0 {
            best = best.min((cnt, cnt2, seed));
            eprintln!(
                "{i:>5} lb {lb} cnt {cnt} {cnt2} extra {} {}   seed {seed}",
                cnt - lb,
                cnt2 - lb
            );
        }
    }
    eprintln!("BEST {best:?}");
    eprintln!(
        "lb {lb} cnt {} {} extra {} {}",
        best.0,
        best.1,
        best.0 - lb,
        best.1 - lb
    );
    eprintln!("lb   density {:.4}", lb as f32 / len as f32);
    eprintln!(
        "best density {:.4} {:.4}",
        best.0 as f32 / len as f32,
        best.1 as f32 / len as f32
    );
    eprintln!(
        "diff density {:.4} {:.4}",
        (best.0 - lb) as f32 / len as f32,
        (best.1 - lb) as f32 / len as f32,
    );
}

// good seeds:
// 4105820918
// 2047504287
// 2090516370
// 650373086
// 303023692 < best?
