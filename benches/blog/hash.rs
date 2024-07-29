pub trait Hasher {
    type Out;
    fn hash(&self, t: &[u8]) -> Self::Out;
    fn hash_kmers(&self, k: usize, t: &[u8]) -> impl Iterator<Item = Self::Out> {
        t.windows(k).map(|kmer| self.hash(kmer))
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FxHash;
impl Hasher for FxHash {
    type Out = u64;
    fn hash(&self, t: &[u8]) -> u64 {
        fxhash::hash64(t)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ExtNtHasher;
impl Hasher for ExtNtHasher {
    type Out = u64;
    fn hash(&self, t: &[u8]) -> u64 {
        nthash::ntf64(t, 0, t.len())
    }
    fn hash_kmers(&self, k: usize, t: &[u8]) -> impl Iterator<Item = Self::Out> {
        nthash::NtHashForwardIterator::new(t, k).unwrap()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct WyHash;
impl Hasher for WyHash {
    type Out = u64;
    fn hash(&self, t: &[u8]) -> u64 {
        wyhash::wyhash(t, 0)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Buffer<H> {
    pub hasher: H,
}
impl<H: Hasher> Hasher for Buffer<H> {
    type Out = H::Out;
    fn hash(&self, t: &[u8]) -> Self::Out {
        self.hasher.hash(t)
    }
    fn hash_kmers(&self, k: usize, t: &[u8]) -> impl Iterator<Item = Self::Out> {
        // Collect the underlying iterator into a Vec before iterating over it.
        self.hasher.hash_kmers(k, t).collect::<Vec<_>>().into_iter()
    }
}
