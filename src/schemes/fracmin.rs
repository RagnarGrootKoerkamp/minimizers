use super::*;

#[derive(Debug, Clone, Serialize)]
pub struct FracMinO {
    bound: u128,
    o: RandomO,
}

impl ToOrder for FracMin {
    type O = FracMinO;
    fn to_order(&self, _w: usize, _k: usize, _sigma: usize) -> FracMinO {
        let bound = u128::MAX / self.f as u128;
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
