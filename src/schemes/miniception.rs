use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MiniceptionNewP {
    pub k0: usize,
}

#[typetag::serialize]
impl Params for MiniceptionNewP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(MiniceptionNew::new(w, k, self.k0, RandomO))
    }
}

/// Modified version of miniception with a different sort order on kmer.
/// Instead of sorting them by h(s), the are ordered by:
/// - First those for which the minimal k0-mer is a prefix, sorted by h(t).
/// - Then those for which the minimal k0-mer is a suffix, sorted by h(t).
pub struct MiniceptionNew<O: Order> {
    w: usize,
    k: usize,
    l: usize,
    k0: usize,
    w0: usize,
    /// The order used for k-mers.
    o: O,
    /// The order used for k0-mers.
    rand_mini: RandomMinimizer,
}

impl<O: Order> MiniceptionNew<O> {
    pub fn new(w: usize, k: usize, k0: usize, o: O) -> Self {
        assert!(k0 >= k.saturating_sub(w));
        assert!(k0 <= k);
        let w0 = k - k0;
        Self {
            w,
            k,
            l: k + w - 1,
            k0,
            w0,
            o,
            rand_mini: RandomMinimizer::new(k0, k - k0 + 1, RandomO),
        }
    }
}

impl<O: Order> SamplingScheme for MiniceptionNew<O> {
    fn l(&self) -> usize {
        self.l
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let o0 = self.rand_mini.ord();
        lmer.windows(self.k)
            .enumerate()
            .filter_map(|(i, kmer)| {
                let j = self.rand_mini.sample(kmer);
                if j == 0 || j == self.k - self.k0 {
                    Some((i, kmer, j == 0, Order::key(o0, &kmer[j..j + self.k0])))
                } else {
                    None
                }
            })
            .min_by_key(|&(_i, w, b, hk0)| (b, hk0, self.o.key(w)))
            .unwrap()
            .0
    }

    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        // Queue of all k0-mers.
        let mut q0 = MonotoneQueue::new();
        // Queue of filtered k-mers.
        let mut q = MonotoneQueue::new();

        // i: position of lmer
        // j: position of kmer
        // j0: position of k0mer

        let o0 = self.rand_mini.ord();

        // 1: init k0-mers.
        let mut k0mers = text.windows(self.k0).enumerate();
        for (j0, k0mer) in k0mers.by_ref().take(self.w0) {
            q0.push(j0, Order::key(o0, k0mer));
        }

        // 2: init k-mers.
        let mut kmers = text.windows(self.k).enumerate().zip(k0mers);
        for ((j, kmer), (j0, k0mer)) in kmers.by_ref().take(self.w - 1) {
            q0.push(j0, Order::key(o0, k0mer));
            let min_pos = q0.pop(j).unwrap().0;
            if min_pos == j || min_pos == j + self.w0 {
                q.push(
                    j,
                    (
                        min_pos == j,
                        Order::key(o0, &text[min_pos..min_pos + self.k0]),
                        self.o.key(kmer),
                    ),
                );
            }
        }

        // 3: Iterate l-mers.
        kmers
            .enumerate()
            .map(move |(i, ((j, kmer), (j0, k0mer)))| {
                q0.push(j0, Order::key(o0, k0mer));
                let min_pos = q0.pop(j).unwrap().0;
                if min_pos == j || min_pos == j + self.w0 {
                    q.push(
                        j,
                        (
                            min_pos == j,
                            Order::key(o0, &text[min_pos..min_pos + self.k0]),
                            self.o.key(kmer),
                        ),
                    );
                }

                q.pop(i).unwrap().0
            })
            .collect()
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct Miniception<O: ToOrder> {
    pub r: usize,
    pub o: O,
}

impl<O: ToOrder> ToOrder for Miniception<O> {
    type O = MiniceptionO;
    fn to_order(&self, _w: usize, k: usize, sigma: usize) -> Self::O {
        let r = self.r;
        MiniceptionO {
            r,
            k,
            m: M::build_from_order(&self.o, k - r + 1, r, sigma),
        }
    }
}

pub struct MiniceptionO {
    k: usize,
    r: usize,
    m: Box<dyn SamplingScheme>,
}

impl Order for MiniceptionO {
    type T = u8;

    fn key(&self, kmer: &[u8]) -> Self::T {
        let x = self.m.sample(kmer);
        if x == 0 || x == self.k - self.r {
            0
        } else {
            1
        }
    }

    fn keys(&self, text: &[u8], _k: usize) -> impl Iterator<Item = Self::T> {
        self.m
            .stream(text)
            .into_iter()
            .map(move |x| if x == 0 || x == self.k - self.r { 0 } else { 1 })
    }
}
