use std::{
    fs,
    io::Write,
    path::PathBuf,
    sync::{atomic::AtomicUsize, Mutex},
};

use clap::Parser;
use itertools::Itertools;
use minimizers::{de_bruijn_seq::de_bruijn_sequence, *};
use order::RandomO;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::Serialize;

#[derive(Debug, Serialize)]
struct Result<'tp> {
    sigma: usize,
    k: usize,
    w: usize,
    l: usize,
    tp: &'tp dyn Params,
    minimizer_name: &'static str,
    density: f64,
    positions: Vec<f64>,
    dists: Vec<f64>,
    transfer: Vec<Vec<f64>>,
}

#[derive(clap::Subcommand)]
enum Command {
    // Run {
    //     #[clap(short, default_value_t = 40)]
    //     w: usize,
    //     #[arg(short, default_value_t = 10)]
    //     k: usize,
    //     #[clap(subcommand)]
    //     tp: MinimizerType,
    // },
    Eval {
        /// Write to file.
        #[clap(short, long)]
        output: Option<PathBuf>,
        /// Run in statistics mode.
        #[clap(short, long)]
        stats: bool,
        #[clap(short, long)]
        small: bool,
        #[clap(short, long)]
        practical: bool,
        #[clap(short, long)]
        practicaly: bool,
        /// Read sequence from file.
        #[clap(short, long)]
        input: Option<PathBuf>,
    },
}

/// Print expected and measured densities for different w and r.
#[derive(clap::Parser)]
struct Args {
    /// Length of the generated random string.
    #[clap(short, default_value_t = 100000)]
    n: usize,
    /// Alphabet size.
    #[clap(short, long, default_value_t = 4)]
    sigma: usize,
    #[clap(subcommand)]
    command: Command,
}

#[ctor::ctor]
fn init_color_backtrace() {
    color_backtrace::install();
}

fn main() {
    let args = Args::parse();
    let sigma = args.sigma;
    let text = &generate_random_string(args.n, sigma);

    match args.command {
        Command::Eval {
            output,
            stats,
            small,
            practical,
            practicaly,
            input,
        } => {
            use schemes::*;
            let r = match sigma {
                2 => 6,
                4 => 4,
                256 => 1,
                _ => panic!("Unsupported alphabet size"),
            };
            let modmini = ModP {
                r,
                t: 0,
                lr: false,
                params: Box::new(M(RandomO)),
            };
            let miniception = schemes::RM(OpenClosed {
                r,
                open: false,
                closed: true,
                open_by_tmer: false,
                closed_by_tmer: false,
                other_by_tmer: false,
                offset: None,
                modulo: false,
                anti_tmer: false,
                o: RandomO,
                miniception_r: true,
            });
            let mut base_types = vec![
                (&schemes::RM(()) as &dyn Params, "RandomMinimizer"),
                // &BdAnchorP { r: 0 },
                (&miniception as &dyn Params, "Miniception"),
                // &MiniceptionNewP { k0: 0 },
                // &ModSamplingP{ k0: 0 },
                (&modmini as &dyn Params, "ModMinimizer"),
                // &ModMinimizerP
                // NOTE: These Rotmini/AltRotmini assume alphabet size 4.
                // &RotMinimizerP
                // &AltRotMinimizerP
                // &DecyclingMinimizerP
                (&schemes::RM(Decycling { double: true }), "DoubleDecycling"),
                // &OpenSyncmerMinimizerP{
                //     t: 0,
                // },
                // &OpenClosedSyncmerMinimizerP { t: 0 },
                // &FracMinP{ f: 0 },
            ];
            if small {
                base_types.push((&BruteforceP as &dyn Params, "Bruteforce"));
            }

            let ks = if small {
                &[1, 2, 3]
            } else if stats {
                &[4, 8, 16, 32]
            } else if practical {
                &(6..=64).collect::<Vec<usize>>()
            } else if practicaly {
                &(6..=32).collect::<Vec<usize>>()
            } else {
                &[
                    5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26,
                    28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 47, 48, 49, 50, 51, 52, 54, 56, 58, 60,
                    62, 64, 66, 68, 70, 72, 74, 76, 78, 80,
                ][..]
            };
            let ws = if small {
                &[2]
            } else if stats {
                &[8, 16, 32, 64]
            } else if practical {
                &[2, 5, 12, 19, 50]
            } else if practicaly {
                &[5, 12, 19]
            } else {
                &[8, 24][..]
            };
            let k_w_tp = ks
                .iter()
                .cartesian_product(ws)
                .cartesian_product(base_types)
                .collect_vec();
            let results = Mutex::new(vec![]);
            let done = AtomicUsize::new(0);
            let total = k_w_tp.len();
            let text = if small {
                &de_bruijn_sequence(sigma, 8)
            } else {
                text
            };
            let text = if let Some(input) = input {
                &fs::read(input).unwrap()
            } else {
                text
            };
            k_w_tp.par_iter().for_each(|&((&k, &w), tp)| {
                if small && k == 3 && sigma == 3 {
                    return;
                }
                let l = w + k - 1;
                let tps = vec![tp];
                for (tp, name) in tps {
                    let (density, positions, dists, transfer, _counts) =
                        collect_stats(w, k, text, &*tp.build(w, k, sigma));
                    let mut results = results.lock().unwrap();
                    results.push(Some(Result {
                        sigma,
                        k,
                        w,
                        l,
                        tp,
                        minimizer_name: name,
                        density,
                        positions: if stats { positions } else { vec![] },
                        dists: if stats { dists } else { vec![] },
                        transfer: if stats { transfer } else { vec![] },
                    }));
                }
                let done = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                eprint!("{done:>3}/{total:>3}: sigma={sigma} k={k} w={w} l={l} tp={tp:?}\r");
            });
            eprintln!();

            if let Some(output) = output {
                let result_json = serde_json::to_string(
                    &results
                        .into_inner()
                        .unwrap()
                        .into_iter()
                        .filter_map(|x| x)
                        .collect_vec(),
                )
                .unwrap();
                let mut file = std::fs::File::create(output).unwrap();
                file.write_all(result_json.as_bytes()).unwrap();
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::{schemes::*, *};
    use order::RandomO;

    #[test]
    fn minimizers() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20 {
            for w in 1..=20 {
                let m = RandomMinimizer::new(k, w, RandomO);
                let stream = m.stream(&text);
                let stream_naive = m.stream_naive(&text);
                assert_eq!(stream, stream_naive);
            }
        }
    }

    #[test]
    fn mod_minimizers() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20 {
            for w in 1..=20 {
                let l = k + w - 1;
                for t in 1..=l {
                    let m = Mod::new_with_t(w, k, 4, t, &M(RandomO));
                    let stream = m.stream(&text);
                    let stream_naive = m.stream_naive(&text);
                    assert_eq!(stream, stream_naive);
                }
            }
        }
    }

    #[test]
    fn miniception() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20usize {
            for w in 1..=20 {
                for k0 in k.saturating_sub(w).max(1)..=k {
                    let m = M((
                        OpenClosed {
                            closed: true,
                            ..OpenClosed::default(k0)
                        },
                        RandomO,
                    ))
                    .build(w, k, 4);
                    let stream = m.stream(&text);
                    let stream_naive = m.stream_naive(&text);
                    assert_eq!(stream, stream_naive);
                }
            }
        }
    }

    #[test]
    fn miniception_new() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20usize {
            for w in 1..=20 {
                for k0 in k.saturating_sub(w).max(1)..=k {
                    let m = M((
                        OpenClosed {
                            closed: true,
                            closed_by_tmer: true,
                            ..OpenClosed::default(k0)
                        },
                        RandomO,
                    ))
                    .build(w, k, 4);
                    let stream = m.stream(&text);
                    let stream_naive = m.stream_naive(&text);
                    assert_eq!(stream, stream_naive);
                }
            }
        }
    }

    #[test]
    fn syncmer_minimizer() {
        let text = generate_random_string(1000, 4);
        for k in 1..=20 {
            for w in 1..=20 {
                for t in 1..k {
                    eprintln!("k {k} w {w} t {t}");
                    let m = M((
                        OpenClosed {
                            open: true,
                            closed: true,
                            open_by_tmer: true,
                            closed_by_tmer: false,
                            ..OpenClosed::default(t)
                        },
                        RandomO,
                    ))
                    .build(w, k, 4);
                    let stream = m.stream(&text);
                    let stream_naive = m.stream_naive(&text);
                    assert_eq!(stream, stream_naive);
                }
            }
        }
    }
}
