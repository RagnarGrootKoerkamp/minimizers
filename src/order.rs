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

/// A directed order: each kmer is mapped to its priority value.
/// Additionally, each kmer is mapped to a direction, determining whether the
/// leftmost or rightmost kmer of this priority is taken in case of ties.
pub trait DirectedOrder {
    fn key(&self, kmer: &[u8]) -> (usize, Direction);
}

/// Maps a kmer to its priority value. Lower is higher priority.
pub trait Order {
    fn key(&self, kmer: &[u8]) -> usize;
}

/// Every order implies a basic directed order.
impl<T: Order> DirectedOrder for T {
    fn key(&self, kmer: &[u8]) -> (usize, Direction) {
        (self.key(kmer), Direction::Leftmost)
    }
}

/// A random order than hashes each kmer using `fxhash64`.
pub struct RandomOrder;

impl Order for RandomOrder {
    fn key(&self, kmer: &[u8]) -> usize {
        fxhash::hash64(kmer) as usize
    }
}

/// An order given by explicit enumeration of all kmers.
#[derive(Clone)]
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

/// Returns the integer value of a kmer.
pub fn pack(kmer: &[u8], sigma: usize) -> usize {
    let mut v = 0;
    for c in kmer {
        v = sigma * v + *c as usize;
    }
    v
}

impl Order for ExplicitOrder {
    fn key(&self, kmer: &[u8]) -> usize {
        assert_eq!(kmer.len(), self.k);
        // Find index of kmer.
        self.idx[pack(kmer, self.sigma)]
    }
}

/// A directed order given by explicit enumeration of all kmers.
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
