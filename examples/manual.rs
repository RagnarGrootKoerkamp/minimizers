#![allow(unused)]
use itertools::Itertools;
use minimizers::{de_bruijn_seq::cyclic_exact_density_string, ExplicitLocalScheme, SamplingScheme};

fn rev(n: usize, k: usize) -> usize {
    let mut n = n;
    let mut r = 0;
    for _ in 0..k {
        r = r * 2 + n % 2;
        n /= 2;
    }
    r
}

fn equivalent(v: &Vec<u8>, k: usize, sigma: usize) -> Vec<Vec<u8>> {
    assert_eq!(sigma, 2);
    let mut eqv = vec![v.clone()];

    // flip 0/1
    {
        let mut rev = v.clone();
        rev.reverse();
        eqv.push(rev);
    }

    // flip left/right
    let mut eqv2 = vec![];
    for v in &eqv {
        let mut v2 = v
            .iter()
            .enumerate()
            .map(|(i, v)| (rev(i, k), k as u8 - 1 - *v))
            .collect_vec();
        v2.sort();
        eqv2.push(v2.iter().map(|(_, v)| *v).collect_vec());
    }
    eqv.extend(eqv2);
    eqv
}

fn contexts(f: &impl SamplingScheme) -> Vec<Vec<u8>> {
    let mut ctxs = vec![];
    let (text, len) = cyclic_exact_density_string(f.k(), f.w(), 2, false);
    let ctx = 2 * f.w() + f.k() - 2;
    let text = &text[..len + ctx - 1];
    let mut seen = std::collections::HashSet::new();
    for (i, pos) in f
        .stream(text)
        .into_iter()
        .skip(f.w() - 1)
        .enumerate()
        .take(len)
    {
        assert!(
            pos - i < ctx - f.k() + 1,
            "{} < {}",
            pos - i,
            ctx - f.k() + 1
        );
        assert!(
            ctx - f.w() - f.k() + 1 <= pos - i,
            "{} <= {}",
            ctx - f.w() - f.k() + 1,
            pos - i
        );
        if seen.contains(&pos) {
            continue;
        }
        seen.insert(pos);

        let mut ctx = &text[i..i + ctx];
        ctxs.push(ctx.to_vec());
        eprintln!("{:>3} / {:>3}: {:?} {}", i, seen.len(), ctx, pos - i);
    }
    ctxs.sort();
    ctxs
}

fn main() {
    let local51 = ExplicitLocalScheme {
        sigma: 2,
        w: 5,
        k: 1,
        map: vec![
            3, 2, 3, 1, 2, 2, 3, 0, 1, 1, 3, 1, 2, 2, 3, 4, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1, 1, 2,
            2, 3, 3,
        ],
    };
    let forward51 = ExplicitLocalScheme {
        sigma: 2,
        w: 5,
        k: 1,
        map: vec![
            3, 2, 3, 1, 2, 2, 3, 0, 1, 1, 1, 1, 2, 2, 3, 4, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1, 1, 2,
            2, 3, 3,
        ],
    };

    let local42 = ExplicitLocalScheme {
        sigma: 2,
        w: 4,
        k: 2,
        map: vec![
            3, 2, 3, 1, 2, 2, 3, 0, 1, 1, 3, 1, 2, 2, 3, 3, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1, 1, 2,
            2, 3, 2,
        ],
    };
    let forward42 = ExplicitLocalScheme {
        sigma: 2,
        w: 4,
        k: 2,
        map: vec![
            3, 2, 3, 1, 2, 2, 3, 0, 1, 1, 1, 1, 2, 2, 3, 3, 0, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1, 1, 2,
            2, 3, 2,
        ],
    };
    let schemes = [&local51, &forward51, &local42, &forward42];

    for f in schemes {
        run(f);
    }

    // let cl = contexts(&local42);
    // let cf = contexts(&forward42);
    let cl = contexts(&local51);
    let cf = contexts(&forward51);
    eprintln!("Only local");
    for c in &cl {
        if !cf.contains(c) {
            eprintln!("{:?}", c);
        }
    }
    eprintln!("Only forward");
    for c in &cf {
        if !cl.contains(c) {
            eprintln!("{:?}", c);
        }
    }
}

fn run(f: &ExplicitLocalScheme) {
    let (text, len) = cyclic_exact_density_string(f.k, f.w, f.sigma, false);
    let cnt = f.cyclic_text_density(&text, len);
    let fwd = f.is_forward(&text);

    eprintln!(
        "{:?}: {:3}/{:3} = {:5.3}, fwd = {}",
        f,
        cnt,
        len,
        cnt as f64 / len as f64,
        fwd
    );
}
