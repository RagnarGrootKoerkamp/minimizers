#![feature(portable_simd)]
#![allow(dead_code)]
mod blog;
use blog::*;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};

criterion_group!(
    name = group;
    config = Criterion::default()
        // Make sure that benchmarks are fast.
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_millis(2000))
        .sample_size(10);
    targets = initial_runtime_comparison,
        blog::counting::count_comparisons_bench,
        optimized, ext_nthash, buffered, local_nthash
);
criterion_main!(group);

fn initial_runtime_comparison(c: &mut Criterion) {
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

fn optimized(c: &mut Criterion) {
    // Create a random string of length 1Mbp.
    let text = &(0..1000000)
        .map(|_| b"ACGT"[rand::random::<u8>() as usize % 4])
        .collect::<Vec<_>>();

    let hasher = FxHash;

    let w = 11;
    let k = 21;

    #[rustfmt::skip]
    let minimizers: &[(&str, &dyn Minimizer, bool)] = &[
        ("buffered_opt", &SlidingWindowMinimizer { w, k, alg: BufferedOpt, hasher }, true),
        ("rescan_opt", &SlidingWindowMinimizer { w, k, alg: RescanOpt, hasher }, true),
        ("split_opt", &SlidingWindowMinimizer { w, k, alg: SplitOpt, hasher }, true),
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

fn ext_nthash(c: &mut Criterion) {
    // Create a random string of length 1Mbp.
    let text = &(0..1000000)
        .map(|_| b"ACGT"[rand::random::<u8>() as usize % 4])
        .collect::<Vec<_>>();

    let w = 11;
    let k = 21;

    let hasher = ExtNtHash;

    #[rustfmt::skip]
    let minimizers: &[(&str, &dyn Minimizer, bool)] = &[
        ("buffered_nt", &SlidingWindowMinimizer { w, k, alg: BufferedOpt, hasher }, true),
        ("queue_nt", &SlidingWindowMinimizer { w, k, alg: Queue, hasher }, true),
        ("jumping_nt", &JumpingMinimizer { w, k, hasher }, false),
        ("rescan_nt", &SlidingWindowMinimizer { w, k, alg: RescanOpt, hasher }, true),
        ("split_nt", &SlidingWindowMinimizer { w, k, alg: SplitOpt, hasher }, true),
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

fn buffered(c: &mut Criterion) {
    // Create a random string of length 1Mbp.
    let text = &(0..1000000)
        .map(|_| b"ACGT"[rand::random::<u8>() as usize % 4])
        .collect::<Vec<_>>();

    let w = 11;
    let k = 21;

    let fxhasher = Buffer { hasher: FxHash };
    let nthasher = Buffer { hasher: ExtNtHash };

    #[rustfmt::skip]
    let minimizers: &[(&str, &dyn Minimizer, bool)] = &[
        ("buffered_buf", &SlidingWindowMinimizer { w, k, alg: BufferedOpt, hasher: fxhasher }, true),
        ("queue_buf", &SlidingWindowMinimizer { w, k, alg: Queue, hasher: fxhasher }, true),
        ("jumping_buf", &JumpingMinimizer { w, k, hasher: fxhasher }, false),
        ("rescan_buf", &SlidingWindowMinimizer { w, k, alg: RescanOpt, hasher: fxhasher }, true),
        ("split_buf", &SlidingWindowMinimizer { w, k, alg: SplitOpt, hasher: fxhasher }, true),
        ("buffered_nt_buf", &SlidingWindowMinimizer { w, k, alg: BufferedOpt, hasher: nthasher }, true),
        ("queue_nt_buf", &SlidingWindowMinimizer { w, k, alg: Queue, hasher: nthasher }, true),
        ("jumping_nt_buf", &JumpingMinimizer { w, k, hasher: nthasher }, false),
        ("rescan_nt_buf", &SlidingWindowMinimizer { w, k, alg: RescanOpt, hasher: nthasher }, true),
        ("split_nt_buf", &SlidingWindowMinimizer { w, k, alg: SplitOpt, hasher: nthasher }, true),
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

fn local_nthash(c: &mut Criterion) {
    // Create a random string of length 1Mbp.
    let text = &(0..1000000)
        .map(|_| b"ACGT"[rand::random::<u8>() as usize % 4])
        .collect::<Vec<_>>();

    let k = 21;

    let mut g = c.benchmark_group("g");
    g.bench_function("ext_nthash", |b| {
        b.iter(|| ExtNtHash.hash_kmers(k, text).collect::<Vec<_>>());
    });
    g.bench_function("nthash", |b| {
        b.iter(|| NtHash.hash_kmers(k, text).collect::<Vec<_>>());
    });
    g.bench_with_input("nthash_buf", text, |b, text| {
        b.iter(|| Buffer { hasher: NtHash }.hash_kmers(k, text));
    });
    g.bench_with_input("nthash_bufopt", text, |b, text| {
        b.iter(|| BufferOpt { hasher: NtHash }.hash_kmers(k, text));
    });
    g.bench_with_input("nthash_bufdouble", text, |b, text| {
        b.iter(|| BufferDouble { hasher: NtHash }.hash_kmers(k, text));
    });
    g.bench_with_input("nthash_bufpar1", text, |b, text| {
        b.iter(|| BufferPar { hasher: NtHashPar::<1> }.hash_kmers(k, text));
    });
    g.bench_with_input("nthash_bufpar2", text, |b, text| {
        b.iter(|| BufferPar { hasher: NtHashPar::<2> }.hash_kmers(k, text));
    });
    g.bench_with_input("nthash_bufpar3", text, |b, text| {
        b.iter(|| BufferPar { hasher: NtHashPar::<3> }.hash_kmers(k, text));
    });
    g.bench_with_input("nthash_bufpar4", text, |b, text| {
        b.iter(|| BufferPar { hasher: NtHashPar::<4> }.hash_kmers(k, text));
    });

    let packed_text = &(0..1000000 / 4)
        .map(|_| rand::random::<u8>())
        .collect::<Vec<_>>();
    g.bench_with_input("nthash_bufsimd", packed_text, |b, packed_text| {
        b.iter(|| BufferPar { hasher: NtHashSimd }.hash_kmers(k, packed_text));
    });
}
