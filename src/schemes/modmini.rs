use super::*;

#[derive(Debug, Serialize)]
pub struct ModP {
    pub r: usize,
    pub lr: bool,
    pub t: usize,
    pub params: Box<dyn Params>,
}

#[typetag::serialize]
impl Params for ModP {
    fn build(&self, w: usize, k: usize, sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(match (self.lr, self.t) {
            (false, 0) => Mod::new(w, k, sigma, self.r, &*self.params),
            (true, _) => Mod::new_lr(w, k, sigma, self.r, &*self.params),
            (_, t) => Mod::new_with_t(w, k, sigma, t, &*self.params),
        })
    }
}

pub struct Mod {
    scheme: Box<dyn SamplingScheme>,
    /// Fast modulo w operations.
    fastmod_w: FM32,
}

impl Mod {
    pub fn new(w: usize, k: usize, sigma: usize, r: usize, params: &dyn Params) -> Self {
        let t = if k < r { k } else { (k - r) % w + r };
        Self::new_with_t(w, k, sigma, t, params)
    }
    pub fn new_lr(w: usize, k: usize, sigma: usize, r: usize, params: &dyn Params) -> Self {
        let t = if k < w + r { k } else { k - w };
        Self::new_with_t(w, k, sigma, t, params)
    }
    pub fn new_with_t(w: usize, k: usize, sigma: usize, t: usize, params: &dyn Params) -> Self {
        let scheme = params.build(w + k - t, t, sigma);
        Self {
            scheme,
            fastmod_w: FM32::new(w),
        }
    }
}

impl SamplingScheme for Mod {
    fn l(&self) -> usize {
        self.scheme.l()
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        self.fastmod_w.reduce(self.scheme.sample(lmer))
    }

    fn stream(&self, text: &[u8]) -> Vec<usize> {
        self.scheme
            .stream(text)
            .into_iter()
            .enumerate()
            .map(|(i, x)| i + self.fastmod_w.reduce(x - i))
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
