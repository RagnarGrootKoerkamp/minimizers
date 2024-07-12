#![feature(duration_constants)]
use itertools::Itertools;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use minimizers::*;
// use order::RandomOrder;

fn bench(c: &mut Criterion) {
    // make a  group
    let mut g = c.benchmark_group("rand");

    let string = generate_random_dna(1000000);
    let k = 30;
    let w = 20;
    // let o = RandomOrder;
    // let mini = RandomMinimizer::new(k, w, o);
    // c.bench_function("simple", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerRehash { k, w, o };
    // c.bench_function("rehash", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerRescan { k, w, o };
    // c.bench_function("rescan", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerRescanNt { k, w };
    // c.bench_function("rescan_nt_0", |b| {
    //     b.iter(|| mini.stream_dedup_0(&string).count());
    // });
    // let mini = MinimizerRescanNt { k, w };
    // c.bench_function("rescan_nt_2", |b| {
    //     b.iter(|| mini.stream_dedup_2(&string).count());
    // });

    // let mini = MinimizerRescanNt::<Tuple>::new(k, w);
    // g.bench_function("rescan_nt_1/tup", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerRescanNt::<Pack>::new(k, w);
    // g.bench_function("rescan_nt_1/pack", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });

    // let mini = MinimizerStacks::<Tuple>::new(k, w);
    // g.bench_function("stacks/tup", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerStacks::<Pack>::new(k, w);
    // g.bench_function("stacks/pack", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerStacksBuf::<Tuple>::new(k, w);
    // g.bench_function("buf/tup", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerStacksBuf::<Pack>::new(k, w);
    // g.bench_function("buf/pack", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });

    let mini = MinimizerStacks::<Pack>::new(k, w);
    // g.bench_function("stacks/pack/dedup", |b| {
    //     b.iter(|| criterion::black_box(mini.stream_dedup_1(&string).collect_vec()));
    // });
    g.bench_function("stacks/pack/all", |b| {
        b.iter(|| criterion::black_box(mini.stream(&string).collect_vec()));
    });

    let mini = MinimizerStacksSimd::new(k, w);
    // g.bench_function("stacks/simd/dedup", |b| {
    //     b.iter(|| criterion::black_box(mini.stream_dedup_1(&string)));
    // });
    g.bench_function("stacks/simd/all", |b| {
        b.iter(|| criterion::black_box(mini.stream(&string)));
    });
}

criterion_group!(
    name = group;
    config = Criterion::default()
        .measurement_time(Duration::SECOND)
        .warm_up_time(Duration::MILLISECOND * 100)
        .sample_size(10);
    targets = bench
);

criterion_main!(group);
