use super::*;

impl ToOrder for RotMinimizer {
    type O = RotMinimizerO;

    fn to_order(&self, w: usize, k: usize, sigma: usize) -> Self::O {
        assert!(k % w == 0);
        RotMinimizerO { w, sigma }
    }
}

impl Order for RotMinimizerO {
    type T = u8;

    fn key(&self, kmer: &[u8]) -> Self::T {
        let psi0 = kmer.iter().step_by(self.w).map(|&x| x as u64).sum::<u64>();
        for j in 1..self.w {
            let psij = kmer
                .iter()
                .skip(j)
                .step_by(self.w)
                .map(|&x| x as u64)
                .sum::<u64>();
            if !(psij <= psi0 + self.sigma as u64) {
                return 1;
            }
        }
        0
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
