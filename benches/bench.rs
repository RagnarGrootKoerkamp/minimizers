#![allow(unused_imports)]
mod randmini;

use criterion::{criterion_group, criterion_main, Criterion};
use itertools::Itertools;
use randmini::{
    baseline::{V0NaiveLex, V1NaiveFx, V2NaiveWy},
    daniel::ExtDaniel,
    queue::{V3AQueue, V3BQueue},
    rescan::V4ARescan,
    Minimizer,
};
use std::time::Duration;

/// Benchmark some functions.
fn bench(c: &mut Criterion) {
    let text = &(0..1000000)
        .map(|_| b"ACGT"[rand::random::<u8>() as usize % 4])
        .collect::<Vec<_>>();
    let k = 21;
    let w = 11;

    let mut g = c.benchmark_group("g");
    // g.bench_function("0_naive_lex", |b| {
    //     let m = V0NaiveLex { w, k };
    //     b.iter(|| m.minimizers(text));
    // });
    // g.bench_function("1_naive_fx", |b| {
    //     let m = V1NaiveFx { w, k };
    //     b.iter(|| m.minimizers(text));
    // });
    // g.bench_function("2_naive_wy", |b| {
    //     let m = V2NaiveWy { w, k };
    //     b.iter(|| m.minimizers(text));
    // });

    g.bench_function("ext_minimizer_iter", |b| {
        b.iter(|| {
            minimizer_iter::MinimizerBuilder::<u64>::new()
                .minimizer_size(k)
                .width(w as u16)
                .iter_pos(text)
                .collect_vec()
        });
    });
    g.bench_function("ext_daniel", |b| {
        let m = ExtDaniel { w, k };
        b.iter(|| m.minimizers(text));
    });
    g.bench_function("3a_queue", |b| {
        let m = V3AQueue { w, k };
        b.iter(|| m.minimizers(text));
    });
    g.bench_function("3b_queue_inline", |b| {
        let m = V3BQueue { w, k };
        b.iter(|| m.minimizers(text));
    });
    g.bench_function("4a_rescan", |b| {
        let m = V4ARescan { w, k };
        b.iter(|| m.minimizers(text));
    });
    g.bench_function("4b_rescan_unchecked", |b| {
        let m = V4ARescan { w, k };
        b.iter(|| m.minimizers(text));
    });
}

criterion_group!(
    name = group;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_millis(2000))
        .sample_size(10);
    targets = bench
);

criterion_main!(group);
