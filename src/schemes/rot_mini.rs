use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RotMinimizerP;

#[typetag::serialize]
impl Params for RotMinimizerP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(RotMinimizer::new(k, w, RandomO))
    }
}

/// Asymptotic (in k) optimal minimizers:
/// - Assume k=x*w; sum all i mod w positions.
/// - Take minimizer if first coordinate is max or at most sigma-1 away from max.
/// - If multiple, take random.
///
/// Assumes alphabet size sigma=4.
/// Also requires the input to be encoded as values 0..sigma, not as arbitrary ascii codepoints.
///
/// On the overshoot:
/// Instead of the original sigma, we can use sigma-1, since the max difference
/// between two characters is sigma-1.  In fact, I conjecture that
/// floor((sigma+1)/2) overshoot is sufficient. I have a sketch of a proof for
/// this.
pub struct RotMinimizer<O: Order> {
    k: usize,
    w: usize,
    o: O,
}
impl<O: Order> RotMinimizer<O> {
    pub fn new(k: usize, w: usize, o: O) -> Self {
        assert!(k % w == 0);
        Self { k, w, o }
    }
}
impl<O: Order> SamplingScheme for RotMinimizer<O> {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let sigma = 4;
        assert!(lmer.iter().all(|c| *c <= sigma));
        lmer.windows(self.k)
            .enumerate()
            .filter(|(_, kmer)| {
                let psi0 = kmer.iter().step_by(self.w).map(|&x| x as u64).sum::<u64>();
                for j in 1..self.w {
                    let psij = kmer
                        .iter()
                        .skip(j)
                        .step_by(self.w)
                        .map(|&x| x as u64)
                        .sum::<u64>();
                    if !(psij <= psi0 + sigma as u64) {
                        return false;
                    }
                }
                true
            })
            .min_by_key(|&(_i, kmer)| self.o.key(kmer))
            .unwrap()
            .0
    }
}

/// Asymptotic (in k) optimal minimizers:
/// - Assume k=x*w; sum all i mod w positions.
/// - Take minimizer if first coordinate is max or at most sigma away from max.
/// - If multiple, take random.
///
/// - Ignore for k < w.
/// - Ignore trailing chars.
///
/// Simplified:
/// - For all k-mers compute the rotational sum of k/w w-strided values.
/// - Take k-mer that maximizes the first coordinate.
#[derive(Clone, Serialize, Debug)]
pub struct AltRotMinimizer;

pub struct AltRotMinimizerO {
    w: usize,
}

impl ToOrder for AltRotMinimizer {
    type O = AltRotMinimizerO;

    fn to_order(&self, w: usize, _k: usize, _sigma: usize) -> Self::O {
        AltRotMinimizerO { w }
    }
}

impl Order for AltRotMinimizerO {
    type T = u64;

    fn key(&self, kmer: &[u8]) -> Self::T {
        kmer.iter().step_by(self.w).map(|&x| x as u64).sum::<u64>()
    }
}
