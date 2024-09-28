use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MinimizerP {
    pub ao: bool,
}

#[typetag::serde]
impl Params for MinimizerP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        if !self.ao {
            Box::new(Minimizer::new(k, w, RandomOrder))
        } else {
            Box::new(Minimizer::new(k, w, AntiLex))
        }
    }
}

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
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|&(_i, kmer)| self.o.key(kmer))
            .unwrap()
            .0
    }

    // TODO: Rolling hash using NtHash.
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        let mut q = MonotoneQueue::new();
        let mut kmers = text.windows(self.k).enumerate();
        // Push the first w-1 kmers onto the queue.
        for (j, kmer) in kmers.by_ref().take(self.w - 1) {
            q.push(j, self.o.key(kmer));
        }
        // i: position of lmer
        // j: position of kmer
        kmers
            .enumerate()
            .map(move |(i, (j, kmer))| {
                q.push(j, self.o.key(kmer));
                q.pop(i).unwrap().0
            })
            .collect()
    }
}
