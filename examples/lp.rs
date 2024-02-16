//! Let X Y and Z be three random variables, such that:
//! - X ~ Y ~ Z, ie they have the same distribution
//! - (X,Y) ~ (Y,Z), ie the pairs have the same joint distribution.
//! - X independent of Z.
//! Then, what is max(P(Y<=X and Y<Z))?
//!
//! This corresponds to k=w=2 minimizers. Given a random string ABCD,
//! X is the hash of AB,
//! Y is the hash of BC,
//! Z is the hash of CD.
//!
//! X, Y, and Z follow the same distribution, as do (X,Y) and (Y,Z).
//! Also X and Z are independent of each other, since AB and CD are disjoint.
//!
//! The probability that the minimizer of BCD is distinct from the minimizer of
//! ABC is the probability that the minimum is either in X or Z, which should be minimized.
//! Here, we maximise the complement of that, ie the probability that Y<=X and Y<Z.
//!
//! For k=w=2 and infinite alphabet:
//! - Random minimizers have density 1/2.
//! - Mod-minimizers have density 5/8.
//! So the question is whether better than this (i.e. P(Y=min) > 3/8) is possible.
//!
//! As it turns out, it is, and in the limit for continuous variables, P(Y=min)=0.5 is possible!
//!
//! This means that:
//! - Either k=w=2 minimizers can achieve density 0.5 (which I think is unlikely),
//! - or the random-variable model lost some assumptions, and the output
//!   distribution here is not actually achievable by a minimizer scheme.

use good_lp::*;

fn main() {
    for t in 1..=5 {
        maximise(t);
    }
}

fn maximise(t: usize) {
    let mut problem = ProblemVariables::new();
    let vars: Vec<Vec<Vec<_>>> = (0..t)
        .map(|_| {
            (0..t)
                .map(|_| {
                    (0..t)
                        .map(|_| problem.add(variable().bounds(0..1)))
                        .collect()
                })
                .collect()
        })
        .collect();

    let mut objective = Expression::default();
    for x in 0..t {
        for y in 0..t {
            for z in 0..t {
                if y <= x && y < z {
                    objective += vars[x][y][z];
                }
            }
        }
    }

    let mut problem = problem.maximise(objective.clone()).using(default_solver);

    // X and Z uniform and independent
    for x in 0..t {
        for z in 0..t {
            let expr: Expression = (0..t).map(|y| vars[x][y][z]).sum();
            problem = problem.with(expr.eq(1. / (t * t) as f64));
        }
    }

    // Y uniform
    for y in 0..t {
        let expr: Expression = (0..t)
            .map(|x| (0..t).map(|z| vars[x][y][z]).sum::<Expression>())
            .sum();
        problem = problem.with(expr.eq(1. / t as f64));
    }

    // (X, Y) ~ (Y, Z)
    for a in 0..t {
        for b in 0..t {
            let e1: Expression = (0..t).map(|z| vars[a][b][z]).sum();
            let e2: Expression = (0..t).map(|x| vars[x][a][b]).sum();
            problem = problem.with(e1.eq(e2));
        }
    }

    let solution = problem.solve().unwrap();
    let out = solution.eval(&objective);
    eprintln!("{t:>3}: {out}");

    for y in 0..t {
        for z in 0..t {
            for x in 0..t {
                eprint!(
                    "{:>.4} ",
                    0.0f64.max(solution.value(vars[x][y][z])) * (t * t * t) as f64
                );
            }
            eprint!("   ");
        }
        eprintln!();
    }
}
