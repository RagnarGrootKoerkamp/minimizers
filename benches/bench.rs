use criterion::{criterion_group, criterion_main, Criterion};
use random::{Minimizer, V0NaiveLex};
use std::time::Duration;

mod random;
use random::*;

/// Benchmark some functions.
fn bench(c: &mut Criterion) {
    let text = &(0..1000000)
        .map(|_| b"ACGT"[rand::random::<u8>() as usize % 4])
        .collect::<Vec<_>>();
    let k = 21;
    let w = 11;

    let mut g = c.benchmark_group("g");
    g.bench_function("0_naive_lex", |b| {
        let m = V0NaiveLex { w, k };
        b.iter(|| m.minimizers(text));
    });
    g.bench_function("1_naive_fx", |b| {
        let m = V1NaiveFx { w, k };
        b.iter(|| m.minimizers(text));
    });
    g.bench_function("2_naive_wy", |b| {
        let m = V2NaiveWy { w, k };
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
