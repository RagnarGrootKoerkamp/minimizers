use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BdAnchorP {
    pub r: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SusAnchorP {
    pub ao: bool,
}

#[typetag::serde]
impl Params for BdAnchorP {
    fn build(&self, w: usize, _k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(BdAnchor::new(w, self.r))
    }
}

#[typetag::serde]
impl Params for SusAnchorP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        if !self.ao {
            Box::new(SusAnchor::new(w, k, Lex))
        } else {
            Box::new(SusAnchor::new(w, k, AntiLex))
        }
    }
}

pub struct BdAnchor {
    w: usize,
    r: usize,
}

impl BdAnchor {
    pub fn new(w: usize, r: usize) -> Self {
        Self { w, r }
    }
}

impl SamplingScheme for BdAnchor {
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
pub struct SusAnchor<O: Order> {
    w: usize,
    k: usize,
    o: O,
}

impl<O: Order> SusAnchor<O> {
    pub fn new(w: usize, k: usize, o: O) -> Self {
        Self { w, k, o }
    }
}

impl<O: Order> SamplingScheme for SusAnchor<O> {
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
