use super::pack;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Ord)]
pub enum Direction {
    Leftmost,
    Rightmost,
}

impl PartialOrd for Direction {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self != other {
            None
        } else if *self == Direction::Leftmost {
            Some(std::cmp::Ordering::Less)
        } else {
            Some(std::cmp::Ordering::Greater)
        }
    }
}

pub trait DirectedOrder {
    fn key(&self, kmer: &[u8]) -> (usize, Direction);
}

pub trait Order {
    fn key(&self, kmer: &[u8]) -> usize;
}

impl<T: Order> DirectedOrder for T {
    fn key(&self, kmer: &[u8]) -> (usize, Direction) {
        (self.key(kmer), Direction::Leftmost)
    }
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

pub struct ExplicitDirectedOrder {
    /// The kmer length of this order.
    pub k: usize,
    /// The alphabet size.
    pub sigma: usize,
    /// the idx for each kmer.
    pub idx: Vec<(usize, Direction)>,
}

impl DirectedOrder for ExplicitDirectedOrder {
    fn key(&self, kmer: &[u8]) -> (usize, Direction) {
        assert_eq!(kmer.len(), self.k);
        // Find index of kmer.
        self.idx[pack(kmer, self.sigma)]
    }
}
