pub trait Hasher {
    type Out;
    fn hash(&self, t: &[u8]) -> Self::Out;
}

#[derive(Clone, Copy, Debug)]
pub struct FxHash;

impl Hasher for FxHash {
    type Out = u64;
    fn hash(&self, t: &[u8]) -> u64 {
        fxhash::hash64(t)
    }
}
