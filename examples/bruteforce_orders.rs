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
use minimizers::{bruteforce, ExplicitLocalScheme};

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
            let best = bruteforce::bruteforce_minimizer(k, w, sigma);
            eprintln!("  best cnt : {} / {}", best.0 .0, best.0 .1);
            eprintln!("  density  : {}", (best.0 .0) as f32 / (best.0 .1) as f32);
            eprintln!("  best perm: {:?}", best.1.ord().idx);
        }
        Method::BFDirectedMinimizer => {
            let best = bruteforce::bruteforce_directed_minimizer(k, w, sigma);
            eprintln!("  best cnt : {} / {}", best.0 .0, best.0 .1);
            eprintln!("  density  : {}", (best.0 .0) as f32 / (best.0 .1) as f32);
            eprintln!("  best perm: {:?}", best.1.ord().idx);
        }
        Method::BFLocalScheme => {
            let best = bruteforce::bruteforce_local_scheme(k, w, sigma);
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
                ExplicitLocalScheme {
                    k: 2,
                    w: 2,
                    sigma: 2,
                    map: vec![0, 1, 0, 0, 1, 1, 0, 1],
                },
                ExplicitLocalScheme {
                    k: 2,
                    w: 2,
                    sigma: 3,
                    map: vec![
                        1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0,
                        0, 1,
                    ],
                },
                ExplicitLocalScheme {
                    k: 2,
                    w: 2,
                    sigma: 4,
                    map: vec![
                        1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1,
                        1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1,
                        0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1,
                    ],
                },
                ExplicitLocalScheme {
                    k: 2,
                    w: 3,
                    sigma: 2,
                    map: vec![0, 2, 1, 1, 0, 2, 0, 0, 1, 2, 1, 1, 2, 2, 1, 2],
                },
                ExplicitLocalScheme {
                    k: 2,
                    w: 3,
                    sigma: 2,
                    map: vec![0, 2, 1, 1, 0, 0, 0, 0, 1, 2, 1, 1, 2, 2, 0, 1],
                },
                ExplicitLocalScheme {
                    k: 2,
                    w: 4,
                    sigma: 2,
                    map: vec![
                        0, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 2, 2, 1, 1, 1, 1, 2,
                        3, 2, 2, 0, 3, 1, 2,
                    ],
                },
                ExplicitLocalScheme {
                    k: 2,
                    w: 4,
                    sigma: 2,
                    map: vec![
                        3, 2, 3, 1, 2, 2, 3, 3, 1, 1, 3, 1, 2, 2, 3, 2, 0, 0, 0, 0, 0, 2, 3, 0, 1,
                        1, 3, 1, 2, 2, 3, 1,
                    ],
                },
                ExplicitLocalScheme {
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
                ExplicitLocalScheme {
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

mod ilp_scheme {

    //! Finds an optimal local scheme using ILP.
    //!
    //! 1. For each lmer, create w boolean variables, 1 of which must be true.
    //! 2. Create a De Bruijn word of sufficient order, and a boolean variable per kmer.
    //! 3. Each kmer selected by the local scheme must be true.
    //! 4. Minimizer the number of selected kmers.

    use itertools::repeat_n;

    use good_lp::*;
    use itertools::Itertools;

    use minimizers::{
        de_bruijn_seq::cyclic_exact_density_string, order::pack, ExplicitLocalScheme,
    };

    pub fn best_local_scheme(
        k: usize,
        w: usize,
        sigma: usize,
        forward: bool,
    ) -> ((usize, usize), ExplicitLocalScheme) {
        let l = k + w - 1;
        let lmers = repeat_n(0..sigma, l)
            .multi_cartesian_product()
            .collect_vec();

        let mut problem = ProblemVariables::new();

        // var[i][j] = 1 if the i'th lmer selects position j.
        let lmer_vars: Vec<Vec<_>> = (0..lmers.len())
            .map(|_| {
                (0..w)
                    .map(|_| problem.add(variable().integer().bounds(0..1)))
                    .collect()
            })
            .collect();

        // len sigma^order + (k-1)
        let (text, chars) = cyclic_exact_density_string(k, w, sigma, forward);
        // eprintln!("#text: {text:?}");
        eprintln!("dBG len: {chars}");
        // double the text.

        let selected_vars = (0..chars)
            .map(|_| problem.add(variable().integer().bounds(0..1)))
            .collect_vec();
        eprintln!("#selected_vars: {}", selected_vars.len());

        let mut objective = Expression::default();
        for selected in &selected_vars {
            objective += selected;
        }

        let mut problem = problem.minimise(objective.clone()).using(default_solver);

        let num_vars = lmer_vars.len() + selected_vars.len();
        eprintln!("#vars: {num_vars}");
        let mut num_constraints = 0;

        // Select 1 per lmer.
        for lmer_var in &lmer_vars {
            problem.add_constraint(lmer_var.iter().sum::<Expression>().eq(1));
            num_constraints += 1;
        }

        // Set selected kmers.
        for (i, lmer) in text.windows(l).take(chars).enumerate() {
            let lmer_id = pack(lmer, sigma);
            for (j, lmer_var) in lmer_vars[lmer_id].iter().enumerate() {
                let selected_var = &selected_vars[(i + j) % chars];
                // eprintln!("{:?} <= {:?}", lmer_var, selected_var);
                problem.add_constraint(lmer_var.into_expression().leq(selected_var));
                num_constraints += 1;
            }
        }

        if forward {
            // let x0..x(w-1) and y0..y(w-1) be the variables for consecutive lmers.
            //
            // x0 x1 x2
            //    y0 y1 y2
            // y0 <= x0 + x1, i.e. if x0=0 and x1=0, then y0=0.
            //
            // Thus:
            // yi <= x0 + x1 + ... + x(i+1) for all 0<=i<w-2. (i=w-2 is redundant.)
            for (lmer1, lmer2) in text.windows(l).tuple_windows().take(chars) {
                let lmer1_id = pack(lmer1, sigma);
                let lmer2_id = pack(lmer2, sigma);
                for i in 0..w - 2 {
                    let sum_xi = lmer_vars[lmer1_id][..=i + 1]
                        .iter()
                        .map(|var| var.into_expression())
                        .sum::<Expression>();
                    let yi = &lmer_vars[lmer2_id][i];
                    problem.add_constraint(yi.into_expression().leq(sum_xi));
                    num_constraints += 1;
                }
            }
        }

        eprintln!("#constraints: {num_constraints}");

        let solution = problem.solve().unwrap();
        let selected = solution.eval(&objective);
        eprintln!("#selected: {selected}");

        let ls = ExplicitLocalScheme {
            k,
            w,
            sigma,
            map: lmer_vars
                .iter()
                .map(|lmer_vars| {
                    lmer_vars
                        .iter()
                        .enumerate()
                        .find(|(_i, var)| solution.value(**var).round() == 1.0)
                        .unwrap()
                        .0 as u8
                })
                .collect(),
        };

        ((selected as usize, chars), ls)
    }
}
