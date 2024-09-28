//! Suppose we have a (w,k) forward sampling scheme.
//! Consider a cycle of length l >= w+k.
//! Let the distances between consecutive sampled positions be d1+d2+...+di = l, with 0 < di <= w.
//! We can consider all partitions D=(di) with their probability p_D of occurring in an l-cycle.
//! From this, we get the full probability distribution over jumped distances.
//! Since l>=w+k, this distribution is the same for all such l.
//! The goal is to minimize the expected jump size.

use clap::Parser;

#[derive(clap::Parser)]
struct Args {
    w: usize,
    k: usize,
    l1: usize,
    l2: usize,
}

fn main() {
    let Args { w, k, l1, l2 } = Args::parse();

    cycle_partition_lp(w, k, l1..l2);
}

use good_lp::*;
use itertools::Itertools;
use std::{iter::zip, ops::Range};

/// Returns all sorted lists of positive numbers that sum to n.
fn partitions(n: usize, max: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![vec![]];
    }
    let mut out = vec![];
    for i in 1..=max.min(n) {
        for mut p in partitions(n - i, i) {
            p.push(i);
            out.push(p);
        }
    }
    out
}

pub fn cycle_partition_lp(w: usize, k: usize, ls: Range<usize>) -> f64 {
    let mut problem = ProblemVariables::new();

    // Add variables for probability of jump distance.
    let p_dist: Vec<_> = (0..=w)
        .map(|_| problem.add(variable().bounds(0..1)))
        .collect();

    let l_w_parts = ls
        .clone()
        .map(|l| {
            // Add weight variable for each partition of l.
            let parts = partitions(l, w);
            eprintln!("Parts: {:?}", parts.len());
            let weights = parts
                .iter()
                .map(|_| problem.add(variable().bounds(0..1)))
                .collect::<Vec<_>>();
            (parts, weights)
        })
        .collect_vec();

    // Maximize expected jump distance.
    let mut objective = Expression::default();
    for d in 1..=w {
        objective += d as f64 * p_dist[d];
    }

    let mut problem = problem.maximise(objective.clone()).using(default_solver);

    // Sum of distance probabilities is 1.
    let s: Expression = p_dist.iter().sum();
    problem.add_constraint(s.eq(1.0));

    // let s: Expression = p_dist[1].into();
    // problem.add_constraint(s.eq(0.0));
    // let s: Expression = p_dist[2].into();
    // problem.add_constraint(s.eq(0.0));

    // For each distance, it is the sum of probabilities of partitions that sum to that distance.

    for (parts, w_parts) in &l_w_parts {
        for d in 1..=w {
            let mut s: Expression = 0.0.into();
            for (part, w_part) in zip(parts, w_parts) {
                let d_count = part.iter().filter(|&&x| x == d).count();
                if d_count > 0 {
                    s += *w_part * d_count as f64;
                }
            }
            problem.add_constraint(s.eq(p_dist[d]));
        }
    }

    let solution = problem.solve().unwrap();
    let avg_dist = solution.eval(&objective);
    let density = 1.0 / avg_dist;
    eprintln!("w{w:>3} k{k:>3} ls{ls:?}: dist={avg_dist} density={density}");
    for d in 0..=w {
        let p = solution.eval(&p_dist[d]);
        eprintln!("dist {d:>2}: {p}");
    }
    for (l, w_parts) in zip(ls, l_w_parts) {
        eprintln!("l={l}");
        for (part, w_part) in zip(w_parts.0, w_parts.1) {
            let p = solution.eval(w_part);
            if p > 0.0 {
                eprintln!(" {p}: {part:?}");
            }
        }
    }
    avg_dist
}
