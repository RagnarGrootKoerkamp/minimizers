use super::*;

#[derive(Debug, Clone, Serialize)]
pub struct Miniception<O: ToOrder> {
    pub r: usize,
    pub by_tmer: bool,
    pub o: O,
}

impl<O: ToOrder> ToOrder for Miniception<O> {
    type O = MiniceptionO<O::O>;
    fn to_order(&self, _w: usize, k: usize, sigma: usize) -> Self::O {
        let r = self.r;
        MiniceptionO {
            r,
            k,
            by_tmer: self.by_tmer,
            m: M::build_from_order(&self.o, k - r + 1, r, sigma),
        }
    }
}

pub struct MiniceptionO<O: Order> {
    k: usize,
    r: usize,
    by_tmer: bool,
    m: Minimizer<O>,
}

impl<O: Order> MiniceptionO<O> {
    fn inner_key(&self, kmer: &[u8], x: usize) -> (u8, O::T) {
        let p = if x == 0 || x == self.k - self.r { 0 } else { 1 };
        let th = if self.by_tmer {
            self.m.ord().key(&kmer[x..x + self.r])
        } else {
            O::T::default()
        };
        (p, th)
    }
}

impl<O: Order> Order for MiniceptionO<O> {
    type T = (u8, O::T);

    fn key(&self, kmer: &[u8]) -> Self::T {
        let x = self.m.sample(kmer);
        self.inner_key(kmer, x)
    }

    fn keys(&self, text: &[u8], k: usize) -> impl Iterator<Item = Self::T> {
        self.m
            .stream(text)
            .into_iter()
            .zip(text.windows(k))
            .map(move |(x, kmer)| self.inner_key(kmer, x))
    }
}
