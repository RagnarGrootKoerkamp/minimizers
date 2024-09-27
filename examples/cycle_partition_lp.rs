use clap::Parser;
use minimizers::cycle_partition_lp::cycle_partition_lp;

#[derive(clap::Parser)]
struct Args {
    w: usize,
    k: usize,
    l1: usize,
    l2: usize,
}

fn main() {
    let Args { w, k, l1, l2 } = Args::parse();

    cycle_partition_lp(w, k, l1..l2);
}
