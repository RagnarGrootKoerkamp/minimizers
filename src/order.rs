use super::pack;

pub trait Order {
    fn key(&self, kmer: &[u8]) -> usize;
}

pub struct RandomOrder;

impl Order for RandomOrder {
    fn key(&self, kmer: &[u8]) -> usize {
        fxhash::hash64(kmer) as usize
    }
}

pub struct ExplicitOrder {
    /// The kmer length of this order.
    pub k: usize,
    /// The alphabet size.
    pub sigma: usize,
    /// the idx for each kmer.
    pub idx: Vec<usize>,
}

impl ExplicitOrder {
    pub fn lexicographic(k: usize, sigma: usize) -> Self {
        let n = sigma.pow(k as u32);
        Self {
            k,
            sigma,
            idx: (0..n).collect(),
        }
    }
}

impl Order for ExplicitOrder {
    fn key(&self, kmer: &[u8]) -> usize {
        assert_eq!(kmer.len(), self.k);
        // Find index of kmer.
        self.idx[pack(kmer, self.sigma)]
    }
}
