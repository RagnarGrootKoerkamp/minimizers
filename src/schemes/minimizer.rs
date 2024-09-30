use super::*;

/// Classic minimizers with respect to some order on kmers `O`.
pub struct Minimizer<O: DirectedOrder> {
    k: usize,
    w: usize,
    l: usize,
    o: O,
}

impl<O: DirectedOrder> Minimizer<O> {
    pub fn new(k: usize, w: usize, o: O) -> Self {
        assert!(k > 0);
        let l = k + w - 1;
        Self { k, w, l, o }
    }
    pub fn ord(&self) -> &O {
        &self.o
    }
}

impl<O: DirectedOrder> SamplingScheme for Minimizer<O> {
    fn w(&self) -> usize {
        self.w
    }
    fn k(&self) -> usize {
        self.k
    }
    fn l(&self) -> usize {
        self.l
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        debug_assert_eq!(lmer.len(), self.l);
        self.o
            .keys(lmer, self.k)
            .enumerate()
            .min_by_key(|&(_i, kmer)| kmer)
            .unwrap()
            .0
    }

    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        let mut q = MonotoneQueue::new();
        let mut ords = self.o.keys(text, self.k).enumerate();
        // Push the first w-1 kmers onto the queue.
        for (j, ord) in ords.by_ref().take(self.w - 1) {
            q.push(j, ord);
        }
        // i: position of lmer
        // j: position of kmer
        ords.enumerate()
            .map(move |(i, (j, ord))| {
                q.push(j, ord);
                q.pop(i).unwrap().0
            })
            .collect()
    }
}
