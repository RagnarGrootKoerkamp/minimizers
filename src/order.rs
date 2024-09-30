#[derive(Clone, Copy, Debug, PartialEq, Eq, Ord)]
pub enum Direction {
    Leftmost,
    Rightmost,
}

impl PartialOrd for Direction {
    #[inline(always)]
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

// TODO: Unify Order and DirectedOrder by returning an associated type O: Ord.

/// A directed order: each kmer is mapped to its priority value.
/// Additionally, each kmer is mapped to a direction, determining whether the
/// leftmost or rightmost kmer of this priority is taken in case of ties.
pub trait DirectedOrder {
    fn key(&self, kmer: &[u8]) -> (impl Ord + Copy, Direction);
    fn keys(&self, text: &[u8], k: usize) -> impl Iterator<Item = (impl Ord + Copy, Direction)> {
        text.windows(k).map(|kmer| self.key(kmer))
    }
}

/// Maps a kmer to its priority value. Lower is higher priority.
pub trait Order {
    fn key(&self, kmer: &[u8]) -> usize;
    fn keys(&self, text: &[u8], k: usize) -> impl Iterator<Item = impl Ord + Copy> {
        text.windows(k).map(|kmer| self.key(kmer))
    }
}

/// Every order implies a basic directed order.
impl<T: Order> DirectedOrder for T {
    #[inline(always)]
    fn key(&self, kmer: &[u8]) -> (impl Ord + Copy, Direction) {
        (self.key(kmer), Direction::Leftmost)
    }
}

/// A random order than hashes each kmer using `wyhash64`.
#[derive(Clone, Copy)]
pub struct RandomOrder;

impl Order for RandomOrder {
    #[inline(always)]
    fn key(&self, kmer: &[u8]) -> usize {
        wyhash::wyhash(kmer, 3141592) as usize
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
    #[inline(always)]
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
#[inline(always)]
pub fn pack(kmer: &[u8], sigma: usize) -> usize {
    let mut v = 0;
    for c in kmer {
        v = sigma * v + *c as usize;
    }
    v
}

impl Order for ExplicitOrder {
    #[inline(always)]
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
    #[inline(always)]
    fn key(&self, kmer: &[u8]) -> (impl Ord + Copy, Direction) {
        assert_eq!(kmer.len(), self.k);
        // Find index of kmer.
        self.idx[pack(kmer, self.sigma)]
    }
}

#[derive(Clone, Copy)]
pub struct Lex;

impl Order for Lex {
    fn key(&self, kmer: &[u8]) -> usize {
        let s = 8usize.saturating_sub(kmer.len());
        let mut prefix = [0xff; 8];
        for i in 0..8 - s {
            prefix[7 - i] = kmer[i];
        }
        usize::from_ne_bytes(prefix)
    }
}

#[derive(Clone, Copy)]
pub struct AntiLex;

impl Order for AntiLex {
    fn key(&self, kmer: &[u8]) -> usize {
        let s = 8usize.saturating_sub(kmer.len());
        let mut prefix = [0xff; 8];
        for i in 0..8 - s {
            prefix[7 - i] = kmer[i];
        }
        prefix[7] ^= 0xff;
        usize::from_ne_bytes(prefix)
    }
}
