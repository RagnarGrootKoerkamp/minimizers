use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpenSyncmerMinimizerP {
    pub t: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpenClosedSyncmerMinimizerP {
    pub t: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OcModMinimizerP {
    pub t: usize,
    pub offset: usize,
    pub use_closed: bool,
    pub prefer_prefix: bool,
    pub open_tmer: bool,
    pub closed_tmer: bool,
    pub other_tmer: bool,
    pub ao: bool,
    pub aot: bool,
}

#[typetag::serde]
impl Params for OpenSyncmerMinimizerP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(OpenSyncmer::new(k, w, self.t, true, false))
    }
}

#[typetag::serde]
impl Params for OpenClosedSyncmerMinimizerP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(OpenSyncmer::new(k, w, self.t, true, true))
    }
}

#[typetag::serde]
impl Params for OcModMinimizerP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        if !self.ao {
            if !self.aot {
                Box::new(OcModMinimizer::new(
                    k,
                    w,
                    self.t,
                    self.offset,
                    self.use_closed,
                    self.prefer_prefix,
                    self.open_tmer,
                    self.closed_tmer,
                    self.other_tmer,
                    RandomOrder,
                    RandomOrder,
                ))
            } else {
                Box::new(OcModMinimizer::new(
                    k,
                    w,
                    self.t,
                    self.offset,
                    self.use_closed,
                    self.prefer_prefix,
                    self.open_tmer,
                    self.closed_tmer,
                    self.other_tmer,
                    RandomOrder,
                    AntiLex,
                ))
            }
        } else {
            if !self.aot {
                Box::new(OcModMinimizer::new(
                    k,
                    w,
                    self.t,
                    self.offset,
                    self.use_closed,
                    self.prefer_prefix,
                    self.open_tmer,
                    self.closed_tmer,
                    self.other_tmer,
                    AntiLex,
                    RandomOrder,
                ))
            } else {
                Box::new(OcModMinimizer::new(
                    k,
                    w,
                    self.t,
                    self.offset,
                    self.use_closed,
                    self.prefer_prefix,
                    self.open_tmer,
                    self.closed_tmer,
                    self.other_tmer,
                    AntiLex,
                    AntiLex,
                ))
            }
        }
    }
}

/// 1. Find open-syncmer kmers, where the smallest tmer is in the middle.
/// 2. Take the minimizer, preferring open-syncmers.
pub struct OpenSyncmer {
    k: usize,
    w: usize,
    /// Length of tmers to consider inside each kmer.
    t: usize,
    o: RandomOrder,
    rand_mini: RandomMinimizer,
    tiebreak: bool,
    closed: bool,
}

impl OpenSyncmer {
    pub fn new(k: usize, w: usize, t: usize, tiebreak: bool, closed: bool) -> Self {
        Self {
            k,
            w,
            t,
            o: RandomOrder,
            rand_mini: RandomMinimizer::new(t, k - t + 1, RandomOrder),
            tiebreak,
            closed,
        }
    }

    fn hash_kmer(&self, kmer: &[u8], j: usize) -> (i32, usize) {
        let m = self.k - self.t;
        let half = m / 2;
        // Prefer kmers whose minimal tmer is in the middle.
        let pref = if j == half {
            0
        } else if self.closed && j == 0 {
            1
        } else if self.closed && j == m {
            1
        } else {
            3
        };
        let tiebreak_hash = if pref == 0 && self.tiebreak {
            let min_tmer = &kmer[j..j + self.t];
            Order::key(self.rand_mini.ord(), min_tmer)
        } else {
            Order::key(&self.o, kmer)
        };
        (pref, tiebreak_hash)
    }
}

impl SamplingScheme for OpenSyncmer {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|(_, kmer)| {
                let j = self.rand_mini.sample(kmer);
                self.hash_kmer(kmer, j)
            })
            .unwrap()
            .0
    }

    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        // Queue of t-mers.
        let mut qt = MonotoneQueue::new();
        // Queue of k-mers.
        let mut q = MonotoneQueue::new();

        // i: position of lmer
        // j: position of kmer
        // jt: position of tmer

        let ot = self.rand_mini.ord();

        // 1: init t-mers.
        let mut tmers = text.windows(self.t).enumerate();
        for (jt, tmer) in tmers.by_ref().take(self.k - self.t) {
            qt.push(jt, Order::key(ot, tmer));
        }

        // 2: init k-mers.
        let mut kmers = text.windows(self.k).enumerate().zip(tmers);
        for ((j, kmer), (jt, tmer)) in kmers.by_ref().take(self.w - 1) {
            qt.push(jt, Order::key(ot, tmer));
            q.push(j, self.hash_kmer(kmer, qt.pop(j).unwrap().0 - j));
        }

        // 3: Iterate l-mers.
        kmers
            .enumerate()
            .map(move |(i, ((j, kmer), (jt, tmer)))| {
                qt.push(jt, Order::key(ot, tmer));
                q.push(j, self.hash_kmer(kmer, qt.pop(j).unwrap().0 - j));
                q.pop(i).unwrap().0
            })
            .collect()
    }
}

/// Prefer first positions t mod w.
/// Then prefer the ends.
pub struct OcModMinimizer<O: Order, OT: Order> {
    k: usize,
    w: usize,
    /// Length of tmers to consider inside each kmer.
    t: usize,
    /// Preferred offset
    offset: usize,
    /// If false, only do the open syncmer part.
    use_closed: bool,
    /// Prefer prefix over suffix?
    prefer_prefix: bool,
    /// Hash middle by tmer.
    open_tmer: bool,
    /// Hash ends by tmer.
    closed_tmer: bool,
    /// Hash remaining cases by tmer
    other_tmer: bool,
    o: O,
    ot: Minimizer<OT>,
}

impl<O: Order, OT: Order> OcModMinimizer<O, OT> {
    pub fn new(
        k: usize,
        w: usize,
        t: usize,
        offset: usize,
        use_closed: bool,
        prefer_prefix: bool,
        open_tmer: bool,
        closed_tmer: bool,
        other_tmer: bool,
        o: O,
        ot: OT,
    ) -> Self {
        Self {
            k,
            w,
            t,
            offset,
            use_closed,
            prefer_prefix,
            open_tmer,
            closed_tmer,
            other_tmer,
            o,
            ot: Minimizer::new(t, k - t + 1, ot),
        }
    }

    fn hash_kmer(&self, kmer: &[u8], j: usize) -> (i32, usize) {
        assert!(j <= self.k - self.t);
        // Prefer kmers whose minimal tmer is in the middle.
        let pref = if j % self.w == self.offset % self.w {
            0
        } else if self.use_closed && j == 0 {
            1
        } else if self.use_closed && j == self.k - self.t {
            if self.prefer_prefix {
                2
            } else {
                1
            }
        } else {
            3
        };
        let hash_tmer = if pref == 0 {
            self.open_tmer
        } else if pref == 1 || pref == 2 {
            self.closed_tmer
        } else {
            self.other_tmer
        };
        let tiebreak_hash = if hash_tmer {
            let min_tmer = &kmer[j..j + self.t];
            Order::key(self.ot.ord(), min_tmer)
        } else {
            Order::key(&self.o, kmer)
        };
        (pref, tiebreak_hash)
    }
}

impl<O: Order, OT: Order> SamplingScheme for OcModMinimizer<O, OT> {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|(_, kmer)| {
                let j = self.ot.sample(kmer);
                self.hash_kmer(kmer, j)
            })
            .unwrap()
            .0
    }

    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        // Queue of t-mers.
        let mut qt = MonotoneQueue::new();
        // Queue of k-mers.
        let mut q = MonotoneQueue::new();

        // i: position of lmer
        // j: position of kmer
        // jt: position of tmer

        let ot = self.ot.ord();

        // 1: init t-mers.
        let mut tmers = text.windows(self.t).enumerate();
        for (jt, tmer) in tmers.by_ref().take(self.k - self.t) {
            qt.push(jt, Order::key(ot, tmer));
        }

        // 2: init k-mers.
        let mut kmers = text.windows(self.k).enumerate().zip(tmers);
        for ((j, kmer), (jt, tmer)) in kmers.by_ref().take(self.w - 1) {
            qt.push(jt, Order::key(ot, tmer));
            q.push(j, self.hash_kmer(kmer, qt.pop(j).unwrap().0 - j));
        }

        // 3: Iterate l-mers.
        kmers
            .enumerate()
            .map(move |(i, ((j, kmer), (jt, tmer)))| {
                qt.push(jt, Order::key(ot, tmer));
                q.push(j, self.hash_kmer(kmer, qt.pop(j).unwrap().0 - j));
                q.pop(i).unwrap().0
            })
            .collect()
    }
}
