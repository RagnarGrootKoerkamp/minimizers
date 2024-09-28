use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModSamplingP {
    pub k0: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LrMinimizerP;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModMinimizerP {
    pub r: usize,
    pub aot: bool,
}

#[typetag::serde]
impl Params for ModSamplingP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(ModSampling::new(k, w, self.k0, RandomOrder))
    }
}

#[typetag::serde]
impl Params for LrMinimizerP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(ModSampling::lr_minimizer(k, w, RandomOrder))
    }
}

#[typetag::serde]
impl Params for ModMinimizerP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        if !self.aot {
            Box::new(ModSampling::mod_minimizer(k, w, self.r, RandomOrder))
        } else {
            Box::new(ModSampling::mod_minimizer(k, w, self.r, AntiLex))
        }
    }
}

pub struct ModSampling<O: Order> {
    w: usize,
    l: usize,
    t: usize,
    /// The number of t-mers in an l-mer.
    w_t: usize,
    /// The order on t-mers.
    o: O,
    /// Fast modulo w operations.
    fastmod_w: FM32,
}

impl<O: Order> ModSampling<O> {
    pub fn new(k: usize, w: usize, t: usize, o: O) -> Self {
        assert!(k > 0);
        assert!(t > 0);
        assert!(w > 0);
        let l = k + w - 1;
        let w_t = l - t + 1;
        Self {
            w,
            l,
            t,
            w_t,
            o,
            fastmod_w: FM32::new(w),
        }
    }
}

impl<O: Order> ModSampling<O> {
    pub fn lr_minimizer(k: usize, w: usize, o: O) -> Self {
        assert!(k > w);
        let t = k - w;
        Self::new(k, w, t, o)
    }
    pub fn mod_minimizer(k: usize, w: usize, r: usize, o: O) -> Self {
        let r = r.min(k);
        let t = (k - r) % w + r;
        Self::new(k, w, t, o)
    }
}

impl<O: Order> SamplingScheme for ModSampling<O> {
    fn l(&self) -> usize {
        self.l
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let x = lmer
            .windows(self.t)
            .enumerate()
            .min_by_key(|&(_i, w)| self.o.key(w))
            .unwrap()
            .0;
        x % self.w
    }

    /// NOTE: This is not always a forward scheme.
    #[inline(always)]
    fn stream(&self, text: &[u8]) -> Vec<usize> {
        let mut q = MonotoneQueue::new();
        let mut tmers = text.windows(self.t).enumerate();

        for (j, tmer) in tmers.by_ref().take(self.w_t - 1) {
            q.push(j, self.o.key(tmer));
        }
        // i: position of lmer
        // j: position of tmer
        tmers
            .enumerate()
            .map(move |(i, (j, tmer))| {
                q.push(j, self.o.key(tmer));
                i + self.fastmod_w.reduce(q.pop(i).unwrap().0 - i)
            })
            .collect()
    }
}

/// Fast %d operations.
/// Taken from https://github.com/lemire/fastmod/blob/master/include/fastmod.h
#[derive(Copy, Clone, Debug)]
struct FM32 {
    modulus: u64,
    inv: u64,
}
impl FM32 {
    fn new(modulus: usize) -> Self {
        assert!(modulus <= u32::MAX as usize);
        Self {
            modulus: modulus as u64,
            inv: (u64::MAX / modulus as u64).wrapping_add(1),
        }
    }
    fn reduce(self, h: usize) -> usize {
        // 64-bit multiplication always gives a 128-bit result.
        // All this casting is just needed to extract the high word; it's not slower in practice.
        let lowbits = self.inv.wrapping_mul(h as u64);
        ((lowbits as u128 * self.modulus as u128) >> 64) as usize
    }
}
