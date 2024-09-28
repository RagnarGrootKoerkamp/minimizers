use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FracMinP {
    pub f: usize,
}

#[typetag::serde]
impl Params for FracMinP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(FracMin::new(k, w, self.f))
    }
}

/// Sample 1/f of the kmers at random, then minimizer to get to 1/w.
pub struct FracMin {
    k: usize,
    w: usize,
    bound: usize,
    o: RandomOrder,
    seed: usize,
}

impl FracMin {
    pub fn new(k: usize, w: usize, f: usize) -> Self {
        let bound = usize::MAX / f;
        Self {
            k,
            w,
            bound,
            o: RandomOrder,
            seed: rand::random(),
        }
    }
}

impl SamplingScheme for FracMin {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|(_, kmer)| {
                let pref = if Order::key(&self.o, kmer) <= self.bound {
                    0
                } else {
                    1
                };
                let kmer_hash = wyhash::wyhash(kmer, self.seed as u64);
                (pref, kmer_hash)
            })
            .unwrap()
            .0
    }
}
