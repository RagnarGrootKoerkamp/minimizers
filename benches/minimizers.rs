#![feature(duration_constants)]
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use minimizers::*;
use order::RandomOrder;

fn bench(c: &mut Criterion) {
    let string = generate_random_dna(1000000);
    let k = 10;
    let w = 10;
    let o = RandomOrder;
    let mini = RandomMinimizer::new(k, w, o);
    c.bench_function("random_minimizer_simple", |b| {
        b.iter(|| mini.stream_dedup_0(&string).count());
    });
}

criterion_group!(
    name = group;
    config = Criterion::default()
        .measurement_time(Duration::SECOND)
        .warm_up_time(Duration::MILLISECOND * 100)
        .sample_size(20);
    targets = bench
);

criterion_main!(group);
