use std::f64::consts::PI;

use super::*;
use num::{complex::Complex64, Zero};

impl ToOrder for Decycling {
    type O = DecyclingO;
    fn to_order(&self, k: usize) -> DecyclingO {
        let roots = (0..k)
            .map(|j| Complex64::exp(2. * PI * Complex64::i() / k as f64 * j as f64))
            .collect();
        DecyclingO {
            k,
            double: self.double,
            roots,
        }
    }
}

#[derive(Clone, Debug)]
pub struct DecyclingO {
    k: usize,
    double: bool,
    roots: Vec<Complex64>,
}

impl Order for DecyclingO {
    type T = usize;
    fn key(&self, kmer: &[u8]) -> usize {
        assert_eq!(kmer.len(), self.k);
        let mut x = Complex64::zero();
        for j in 0..self.k {
            x += kmer[j] as f64 * self.roots[j];
        }
        let arg = x.arg();
        if PI - 2. * PI / self.k as f64 <= arg && arg < PI {
            0
        } else if self.double && -2. * PI / self.k as f64 <= arg && arg < 0. {
            1
        } else {
            2
        }
    }
}
