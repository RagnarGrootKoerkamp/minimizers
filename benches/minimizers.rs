#![feature(duration_constants)]
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use minimizers::*;
// use order::RandomOrder;

fn bench(c: &mut Criterion) {
    let string = generate_random_dna(1000000);
    let k = 30;
    let w = 20;
    // let o = RandomOrder;
    // let mini = RandomMinimizer::new(k, w, o);
    // c.bench_function("minimizer_simple", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerRehash { k, w, o };
    // c.bench_function("minimizer_rehash", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerRescan { k, w, o };
    // c.bench_function("minimizer_rescan", |b| {
    //     b.iter(|| mini.stream_dedup_1(&string).count());
    // });
    // let mini = MinimizerRescanNt { k, w };
    // c.bench_function("minimizer_rescan_nt_0", |b| {
    //     b.iter(|| mini.stream_dedup_0(&string).count());
    // });
    // let mini = MinimizerRescanNt { k, w };
    // c.bench_function("minimizer_rescan_nt_2", |b| {
    //     b.iter(|| mini.stream_dedup_2(&string).count());
    // });

    let mini = MinimizerRescanNt::<Tuple>::new(k, w);
    c.bench_function("minimizer_rescan_nt_1/tup", |b| {
        b.iter(|| mini.stream_dedup_1(&string).count());
    });
    let mini = MinimizerRescanNt::<Pack>::new(k, w);
    c.bench_function("minimizer_rescan_nt_1/pack", |b| {
        b.iter(|| mini.stream_dedup_1(&string).count());
    });

    let mini = MinimizerStacks::<Tuple>::new(k, w);
    c.bench_function("minimizer_stacks/tup", |b| {
        b.iter(|| mini.stream_dedup_1(&string).count());
    });
    let mini = MinimizerStacks::<Pack>::new(k, w);
    c.bench_function("minimizer_stacks/pack", |b| {
        b.iter(|| mini.stream_dedup_1(&string).count());
    });
    let mini = MinimizerStacksBuf::<Tuple>::new(k, w);
    c.bench_function("minimizer_buf/tup", |b| {
        b.iter(|| mini.stream_dedup_1(&string).count());
    });
    let mini = MinimizerStacksBuf::<Pack>::new(k, w);
    c.bench_function("minimizer_buf/pack", |b| {
        b.iter(|| mini.stream_dedup_1(&string).count());
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
