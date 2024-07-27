use super::*;

#[derive(Clone, Copy, Ord, PartialOrd, Eq, PartialEq)]
pub struct Elem<V> {
    pub val: V,
    pub pos: usize,
}

pub trait SlidingMin<V> {
    /// Take an iterator over values of type V.
    /// Return an iterator over the minima of windows of size w, and their positions.
    fn sliding_min(w: usize, it: impl Iterator<Item = V>) -> impl Iterator<Item = Elem<V>>;
}

/// A iterator extension trait so we can conveniently call .sliding_min on any iterator.
pub trait SlidingMinExtension<V> {
    fn sliding_min<Alg: SlidingMin<V>>(self, w: usize) -> impl Iterator<Item = Elem<V>>;
}

impl<V, I> SlidingMinExtension<V> for I
where
    I: Iterator<Item = V>,
{
    fn sliding_min<Alg: SlidingMin<V>>(self, w: usize) -> impl Iterator<Item = Elem<V>> {
        Alg::sliding_min(w, self)
    }
}

/// A minimizer implementation based on a sliding window minimum algorithm.
pub struct SlidingMinMinimizer<SlidingMinAlg> {
    pub w: usize,
    pub k: usize,
    /// Rust requires that the Alg type parameter is 'used'.
    pub alg: PhantomData<SlidingMinAlg>,
}

impl<SlidingMinAlg: SlidingMin<u64>> Minimizer for SlidingMinMinimizer<SlidingMinAlg> {
    fn window_minimizers(&self, text: &[u8]) -> Vec<usize> {
        // Iterate over k-mers, hash them, and take sliding window minima.
        text.windows(self.k)
            .map(|kmer| fxhash::hash64(kmer))
            // (See source for the iterator 'extension trait'.)
            .sliding_min::<SlidingMinAlg>(self.w)
            .map(|elem| elem.pos)
            .collect()
    }
}

/// A type alias for the corresponding minimizer scheme.
pub type Buffered = SlidingMinMinimizer<BufferedSlidingMin>;
