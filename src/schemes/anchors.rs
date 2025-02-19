use super::*;

#[typetag::serialize]
impl Params for BdAnchor {
    fn build(&self, w: usize, _k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(BdAnchorS::new(w, self.r))
    }
}

#[typetag::serialize]
impl<O: ToOrder> Params for SusAnchor<O> {
    fn build(&self, w: usize, k: usize, sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(SusAnchorS::new(w, k, self.0.to_order(w, k, sigma)))
    }
}

pub struct BdAnchorS {
    w: usize,
    r: usize,
}

impl BdAnchorS {
    pub fn new(w: usize, r: usize) -> Self {
        Self { w, r }
    }
}

impl SamplingScheme for BdAnchorS {
    fn l(&self) -> usize {
        self.w
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let w = self.w;
        debug_assert_eq!(lmer.len(), w);
        let mut best = 0;
        for i in 1..w.saturating_sub(self.r) {
            for j in 0..w {
                if lmer[(i + j) % w] != lmer[(best + j) % w] {
                    if lmer[(i + j) % w] < lmer[(best + j) % w] {
                        best = i;
                    }
                    break;
                }
            }
        }
        best
    }
}

/// NOTE: O should be Lex or AntiLex order. Random order will not be good.
pub struct SusAnchorS<O: Order> {
    w: usize,
    k: usize,
    o: O,
}

impl<O: Order> SusAnchorS<O> {
    pub fn new(w: usize, k: usize, o: O) -> Self {
        Self { w, k, o }
    }
}

impl<O: Order> SamplingScheme for SusAnchorS<O> {
    fn l(&self) -> usize {
        self.w + self.k - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let mut best = (self.o.key(lmer), 0);
        for i in 1..self.w {
            best = best.min((self.o.key(&lmer[i..]), i));
        }
        best.1
    }
}
