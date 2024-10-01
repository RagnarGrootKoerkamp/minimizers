use super::*;

#[derive(Debug, Clone, Serialize)]
pub struct OpenClosed<O: ToOrder> {
    pub r: usize,
    /// Do open syncmers first?
    pub open: bool,
    /// Do closed syncmers second?
    pub closed: bool,
    /// When true, open syncmers must have this offset instead of in the middle by default.
    pub offset: Option<usize>,
    /// When true, any position offset%w makes a kmer an open syncmer.
    pub modulo: bool,
    // FIXME: tiebreak by reverse tmer.
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
            offset: None,
            modulo: false,
            o: RandomO,
        }
    }
}

impl<O: ToOrder> ToOrder for OpenClosed<O> {
    type O = OpenClosedO<O::O>;
    fn to_order(&self, w: usize, k: usize, sigma: usize) -> Self::O {
        let Self {
            r,
            open,
            closed,
            open_by_tmer,
            closed_by_tmer,
            other_by_tmer,
            offset,
            modulo,
            ..
        }: OpenClosed<O> = *self;
        let mut offset = self.offset.unwrap_or((k - r) / 2);
        if modulo {
            offset %= w;
        }
        OpenClosedO {
            r,
            w,
            k,
            open,
            closed,
            open_by_tmer,
            closed_by_tmer,
            other_by_tmer,
            offset,
            modulo,
            m: M::build_from_order(&self.o, k - r + 1, r, sigma),
        }
    }
}

pub struct OpenClosedO<O: Order> {
    r: usize,
    w: usize,
    k: usize,
    open: bool,
    closed: bool,
    open_by_tmer: bool,
    closed_by_tmer: bool,
    other_by_tmer: bool,
    offset: usize,
    modulo: bool,
    m: Minimizer<O>,
}

impl<O: Order> OpenClosedO<O> {
    fn inner_key(&self, kmer: &[u8], x: usize) -> (u8, O::T) {
        let w0 = self.k - self.r;
        let p;
        let by_tmer;
        // FIXME: Re-use fastreduce.
        let is_open = if self.modulo { x % self.w } else { x } == self.offset;
        let is_closed = x == 0 || x == w0;
        if self.open && is_open {
            p = 0;
            by_tmer = self.open_by_tmer;
        } else if self.closed && is_closed {
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
