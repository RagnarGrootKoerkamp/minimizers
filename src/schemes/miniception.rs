use super::*;

#[derive(Debug, Clone, Serialize)]
pub struct OpenClosed<O: ToOrder> {
    pub r: usize,
    pub open: bool,
    pub closed: bool,
    // TODO: By reverse tmer.
    pub open_by_tmer: bool,
    pub closed_by_tmer: bool,
    pub other_by_tmer: bool,
    pub o: O,
}

impl OpenClosed<RandomO> {
    pub fn default(r: usize) -> Self {
        Self {
            r,
            open: false,
            closed: false,
            open_by_tmer: false,
            closed_by_tmer: false,
            other_by_tmer: false,
            o: RandomO,
        }
    }
}

impl<O: ToOrder> ToOrder for OpenClosed<O> {
    type O = OpenClosedO<O::O>;
    fn to_order(&self, _w: usize, k: usize, sigma: usize) -> Self::O {
        let Self {
            r,
            open,
            closed,
            open_by_tmer,
            closed_by_tmer,
            other_by_tmer,
            ..
        } = *self;
        OpenClosedO {
            r,
            k,
            open,
            closed,
            open_by_tmer,
            closed_by_tmer,
            other_by_tmer,
            m: M::build_from_order(&self.o, k - r + 1, r, sigma),
        }
    }
}

pub struct OpenClosedO<O: Order> {
    k: usize,
    r: usize,
    open: bool,
    closed: bool,
    open_by_tmer: bool,
    closed_by_tmer: bool,
    other_by_tmer: bool,
    m: Minimizer<O>,
}

impl<O: Order> OpenClosedO<O> {
    fn inner_key(&self, kmer: &[u8], x: usize) -> (u8, O::T) {
        let w0 = self.k - self.r;
        let half = w0 / 2;
        let p;
        let by_tmer;
        if self.open && x == half {
            p = 0;
            by_tmer = self.open_by_tmer;
        } else if self.closed && (x == 0 || x == w0) {
            p = 1;
            by_tmer = self.closed_by_tmer;
        } else {
            p = 2;
            by_tmer = self.other_by_tmer;
        }
        let tiebreak = if by_tmer {
            self.m.ord().key(&kmer[x..x + self.r])
        } else {
            O::T::default()
        };
        (p, tiebreak)
    }
}

impl<O: Order> Order for OpenClosedO<O> {
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
