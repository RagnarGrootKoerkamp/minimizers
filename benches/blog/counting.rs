use super::*;
use std::cell::Cell;
use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CountCompare<'a, V> {
    val: V,
    count_cmp: Option<&'a Cell<usize>>,
}

impl<V: Max> Max for CountCompare<'_, V> {
    const MAX: Self = CountCompare { val: V::MAX, count_cmp: None };
}

impl<V: Ord> PartialOrd for CountCompare<'_, V> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<V: Ord> Ord for CountCompare<'_, V> {
    fn cmp(&self, other: &Self) -> Ordering {
        // Do not count comparisons with MAX sentinel values.
        if self.count_cmp.is_some() && other.count_cmp.is_some() {
            self.count_cmp
                .unwrap()
                .set(self.count_cmp.unwrap().get() + 1);
        }
        self.val.cmp(&other.val)
    }
}

#[derive(Clone, Copy)]
struct CountingHash<'a, H: Hasher> {
    count_cmp: &'a Cell<usize>,
    hasher: H,
}

impl<'a, H: Hasher> Hasher for CountingHash<'a, H> {
    type Out = CountCompare<'a, H::Out>;
    fn hash(&self, t: &[u8]) -> Self::Out {
        CountCompare { val: self.hasher.hash(t), count_cmp: Some(self.count_cmp) }
    }
}

/// Count the number of operations for each minimizer algorithm.
/// Invoke using `cargo test --bench bench count_comparisons`.
pub fn count_comparisons(_c: &mut criterion::Criterion) {
    let k = 21;
    let n = 100000000 + k - 1;

    for w in [10, 20] {
        let text = &(0..n)
            .map(|_| b"ACGT"[rand::random::<u8>() as usize % 4])
            .collect::<Vec<_>>();
        let n = n as f32;

        let cnt = Cell::new(0);
        let hasher = CountingHash { count_cmp: &cnt, hasher: FxHash };

        let minimizers: &[(&str, &dyn Minimizer); 5] = &[
            (
                "Buffered",
                &SlidingWindowMinimizer { w, k, alg: BufferedSlidingMin, hasher },
            ),
            (
                "Queue",
                &SlidingWindowMinimizer { w, k, alg: Queue, hasher },
            ),
            (
                "Rescan",
                &SlidingWindowMinimizer { w, k, alg: Rescan, hasher },
            ),
            (
                "Split",
                &SlidingWindowMinimizer { w, k, alg: Split, hasher },
            ),
            (
                "Rescan2",
                &SlidingWindowMinimizer { w, k, alg: Rescan2, hasher },
            ),
        ];

        for (name, minimizer) in minimizers.iter() {
            cnt.set(0);
            minimizer.window_minimizers(text);
            println!("{name:<10}: w={w:>2} {:.3}", cnt.get() as f32 / n);
        }

        let (minimizer, name) = (&JumpingMinimizer { w, k, hasher }, "Jumping");
        cnt.set(0);
        minimizer.minimizer_positions(text);
        println!("{name:<10}: w={w:>2} {:.3}", cnt.get() as f32 / n);
    }
}
