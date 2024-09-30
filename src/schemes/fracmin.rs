use super::*;

#[derive(Debug, Clone, Serialize)]
pub struct FracMinO {
    bound: usize,
    o: RandomO,
}

impl ToOrder for FracMin {
    type O = FracMinO;
    fn to_order(&self, _w: usize, _k: usize, _sigma: usize) -> FracMinO {
        let bound = usize::MAX / self.f;
        FracMinO { bound, o: RandomO }
    }
}

impl Order for FracMinO {
    type T = u8;
    fn key(&self, kmer: &[u8]) -> u8 {
        if Order::key(&self.o, kmer) <= self.bound {
            0
        } else {
            1
        }
    }
}
