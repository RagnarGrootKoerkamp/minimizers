#![allow(refining_impl_trait)]

use serde::Serialize;
use std::fmt::Debug;

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
pub trait DirectedOrder: Sync {
    type T: Ord + Copy;
    fn key(&self, kmer: &[u8]) -> (Self::T, Direction);
    fn keys(&self, text: &[u8], k: usize) -> impl Iterator<Item = (Self::T, Direction)> {
        text.windows(k).map(|kmer| self.key(kmer))
    }
}

/// Maps a kmer to its priority value. Lower is higher priority.
pub trait Order: Sync {
    type T: Ord + Copy;
    fn key(&self, kmer: &[u8]) -> Self::T;
    fn keys(&self, text: &[u8], k: usize) -> impl Iterator<Item = Self::T> {
        text.windows(k).map(|kmer| self.key(kmer))
    }
}

/// Every order implies a basic directed order.
impl<O: Order> DirectedOrder for O {
    type T = O::T;
    #[inline(always)]
    fn key(&self, kmer: &[u8]) -> (Self::T, Direction) {
        (self.key(kmer), Direction::Leftmost)
    }
}

/// A random order than hashes each kmer using `wyhash64`.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct RandomO;

impl Order for RandomO {
    type T = usize;
    #[inline(always)]
    fn key(&self, kmer: &[u8]) -> usize {
        wyhash::wyhash(kmer, 3141592) as usize
    }
}

/// An order given by explicit enumeration of all kmers.
#[derive(Clone, Debug)]
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
    type T = usize;
    #[inline(always)]
    fn key(&self, kmer: &[u8]) -> usize {
        assert_eq!(kmer.len(), self.k);
        // Find index of kmer.
        self.idx[pack(kmer, self.sigma)]
    }
}

/// A directed order given by explicit enumeration of all kmers.
#[derive(Clone, Debug)]
pub struct ExplicitDirectedOrder {
    /// The kmer length of this order.
    pub k: usize,
    /// The alphabet size.
    pub sigma: usize,
    /// the idx for each kmer.
    pub idx: Vec<(usize, Direction)>,
}

impl DirectedOrder for ExplicitDirectedOrder {
    type T = usize;
    #[inline(always)]
    fn key(&self, kmer: &[u8]) -> (usize, Direction) {
        assert_eq!(kmer.len(), self.k);
        // Find index of kmer.
        self.idx[pack(kmer, self.sigma)]
    }
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Lex;

impl Order for Lex {
    type T = usize;
    fn key(&self, kmer: &[u8]) -> usize {
        let s = 8usize.saturating_sub(kmer.len());
        let mut prefix = [0xff; 8];
        for i in 0..8 - s {
            prefix[7 - i] = kmer[i];
        }
        usize::from_ne_bytes(prefix)
    }
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct AntiLex;

impl Order for AntiLex {
    type T = usize;
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

impl<O1: Order, O2: Order> Order for (O1, O2) {
    type T = (O1::T, O2::T);
    fn key(&self, kmer: &[u8]) -> Self::T {
        let (o1, o2) = self;
        (o1.key(kmer), o2.key(kmer))
    }

    fn keys(&self, text: &[u8], k: usize) -> impl Iterator<Item = Self::T> {
        let keys0 = self.0.keys(text, k);
        let keys1 = self.1.keys(text, k);
        keys0.zip(keys1)
    }
}

pub trait ToOrder: Clone + Debug + Sync + Serialize + 'static {
    type O: Order;
    fn to_order(&self, w: usize, k: usize, sigma: usize) -> Self::O;
}

impl<O1: ToOrder, O2: ToOrder> ToOrder for (O1, O2) {
    type O = (O1::O, O2::O);
    fn to_order(&self, w: usize, k: usize, sigma: usize) -> Self::O {
        let (o1, o2) = self;
        (o1.to_order(w, k, sigma), o2.to_order(w, k, sigma))
    }
}

impl ToOrder for RandomO {
    type O = RandomO;
    fn to_order(&self, _w: usize, _k: usize, _sigma: usize) -> Self::O {
        RandomO
    }
}

impl ToOrder for Lex {
    type O = Lex;
    fn to_order(&self, _w: usize, _k: usize, _sigma: usize) -> Self::O {
        Lex
    }
}

impl ToOrder for AntiLex {
    type O = AntiLex;
    fn to_order(&self, _w: usize, _k: usize, _sigma: usize) -> Self::O {
        AntiLex
    }
}
