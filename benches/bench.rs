use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};

/// Benchmark some functions.
fn bench(c: &mut Criterion) {
    let text = &(0..1000000)
        .map(|_| b"ACGT"[rand::random::<u8>() as usize % 4])
        .collect::<Vec<_>>();
    let mut g = c.benchmark_group("g");
    g.bench_function("sum of chars", |b| {
        b.iter(|| text.iter().map(|&c| c as usize).sum::<usize>() as usize);
    });
    g.bench_function("sum of squares", |b| {
        b.iter(|| text.iter().map(|&c| (c as usize).pow(2)).sum::<usize>() as usize);
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
