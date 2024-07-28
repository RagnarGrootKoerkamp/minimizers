#![allow(dead_code)]
mod blog;
use blog::*;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};

pub fn initial_runtime_comparison(c: &mut Criterion) {
    // Create a random string of length 1Mbp.
    let text = &(0..1000000)
        .map(|_| b"ACGT"[rand::random::<u8>() as usize % 4])
        .collect::<Vec<_>>();

    let hasher = FxHash;

    let w = 11;
    let k = 21;

    #[rustfmt::skip]
    let minimizers: &[(&str, &dyn Minimizer, bool)] = &[
        ("naive", &NaiveMinimizer { w, k, hasher }, true),
        ("buffered", &SlidingWindowMinimizer { w, k, alg: Buffered, hasher }, true),
        ("queue", &SlidingWindowMinimizer { w, k, alg: Queue, hasher }, true),
        ("jumping", &JumpingMinimizer { w, k, hasher }, false),
        ("rescan", &SlidingWindowMinimizer { w, k, alg: Rescan, hasher }, true),
        ("split", &SlidingWindowMinimizer { w, k, alg: Split, hasher }, true),
        ("queue_igor", &QueueIgor { w, k }, false),
        ("rescan_daniel", &RescanDaniel { w, k }, true),
    ];

    let mut g = c.benchmark_group("g");
    for &(name, m, window_minimizers) in minimizers {
        g.bench_function(name, |b| {
            if window_minimizers {
                b.iter(|| m.window_minimizers(text));
            } else {
                b.iter(|| m.minimizer_positions(text));
            }
        });
    }
}

criterion_group!(
    name = group;
    config = Criterion::default()
        // Make sure that benchmarks are fast.
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_millis(2000))
        .sample_size(10);
    targets = initial_runtime_comparison, blog::counting::count_comparisons_bench
);
criterion_main!(group);
