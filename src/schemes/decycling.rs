use std::f64::consts::PI;

use super::*;
use num::{complex::Complex64, Zero};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DecyclingP;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DoubleDecyclingP {
    pub ao: bool,
}

#[typetag::serde]
impl Params for DecyclingP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(Decycling::new(k, w, RandomOrder, false))
    }
}

#[typetag::serde]
impl Params for DoubleDecyclingP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        if !self.ao {
            Box::new(Decycling::new(k, w, RandomOrder, true))
        } else {
            Box::new(Decycling::new(k, w, AntiLex, true))
        }
    }
}

pub struct Decycling<O: Order> {
    k: usize,
    w: usize,
    roots: Vec<Complex64>,
    o: O,
    double_decycling: bool,
}

impl<O: Order> Decycling<O> {
    pub fn new(k: usize, w: usize, o: O, double_decycling: bool) -> Self {
        let roots = (0..k)
            .map(|j| Complex64::exp(2. * PI * Complex64::i() / k as f64 * j as f64))
            .collect();
        Self {
            k,
            w,
            roots,
            o,
            double_decycling,
        }
    }
}

impl<O: Order> SamplingScheme for Decycling<O> {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    // TODO: Implement streaming version.
    fn sample(&self, lmer: &[u8]) -> usize {
        lmer.windows(self.k)
            .enumerate()
            .min_by_key(|(_i, kmer)| {
                let mut x = Complex64::zero();
                for j in 0..self.k {
                    x += kmer[j] as f64 * self.roots[j];
                }
                let arg = x.arg();
                let key = self.o.key(kmer);
                if PI - 2. * PI / self.k as f64 <= arg && arg < PI {
                    (0, key)
                } else if self.double_decycling && -2. * PI / self.k as f64 <= arg && arg < 0. {
                    (1, key)
                } else {
                    (2, key)
                }
            })
            .unwrap()
            .0
    }
}
