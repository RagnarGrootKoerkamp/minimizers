use std::cell::Cell;

use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThresholdMinimizerP {
    pub t: usize,
    pub h: usize,

    pub loose: bool,
    pub open: bool,
}

#[typetag::serialize]
impl Params for ThresholdMinimizerP {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(ThresholdMinimizer::new(
            k, w, self.t, self.loose, self.open, self.h,
        ))
    }
}
pub struct ThresholdMinimizer {
    k: usize,
    w: usize,
    /// Length of tmers to consider inside each kmer.
    t: usize,
    h: usize,
    loose: bool,
    open: bool,
    o: RandomO,
    rand_mini: RandomMinimizer,
    threshold: usize,
    threshold2: usize,
    count: Cell<usize>,
    count2: Cell<usize>,
    count3: Cell<usize>,
}

impl ThresholdMinimizer {
    pub fn new(k: usize, w: usize, t: usize, loose: bool, open: bool, h: usize) -> Self {
        // let k2 = k - t + 1;
        let threshold = t * h;
        let threshold2 = (2 * t - 256) * h;
        eprintln!("Thresholds: {} {}", threshold, threshold2);
        Self {
            k,
            w,
            t,
            h,
            loose,
            open,
            o: RandomO,
            // rand_mini: Minimizer::new(t, k - t + 1, RandomOrder),
            rand_mini: RandomMinimizer::new(1, k, RandomO),
            threshold,
            threshold2,
            count: Cell::new(0),
            count2: Cell::new(0),
            count3: Cell::new(0),
        }
    }

    fn hash_kmer(&self, kmer: &[u8]) -> (i32, isize, usize) {
        // Hash all tmers in the kmer.
        // If the first tmer is the only one > threshold, it's a good kmer.
        let mut pref = 0;
        // if Order::key(&self.rand_mini.o, &kmer[0..self.t]) >= self.threshold {
        //     for tmer in kmer[1..].windows(self.t) {
        //         if Order::key(&self.rand_mini.o, tmer) < self.threshold {
        //             pref = 1;
        //             break;
        //         }
        //     }
        // } else {
        //     pref = 1;
        // }
        pref = 100;
        let e = self.h.min(self.k);

        let s0: usize = kmer[0..e].iter().map(|&c| c as usize).sum();
        if s0 >= self.threshold {
            pref -= 1;
        }

        if pref < 100 {
            for hmer in kmer[1..].windows(self.h) {
                if hmer.iter().map(|&c| c as usize).sum::<usize>() < self.threshold {
                    pref -= 1;
                } else {
                    break;
                }
            }
        }

        if false {
            // Right-closed kmers
            // if pref > 100 - self.k as i32 {
            //     pref = 200;
            //     if (kmer[self.k - 1] as usize) >= self.threshold {
            //         pref -= 1;
            //     }

            //     if pref < 200 {
            //         for hmer in kmer[..self.k - 1].windows(self.h).rev() {
            //             if hmer.iter().map(|&c| c as usize).sum::<usize>() < self.threshold {
            //                 pref -= 1;
            //             } else {
            //                 break;
            //             }
            //         }
            //     }

            //     if !self.loose && pref > 200 - self.k as i32 {
            //         pref = 200;
            //     }
            // }
        } else {
            if !self.loose && pref > 100 - self.k as i32 {
                pref = 100;
            }
        }

        // if self.open && pref == 100 {
        //     pref = 200;
        //     let s0: usize = kmer[0..e].iter().map(|&c| c as usize).sum();
        //     if s0 >= self.threshold2 {
        //         pref -= 1;
        //     }

        //     if pref < 200 {
        //         for hmer in kmer[1..].windows(self.h) {
        //             if hmer.iter().map(|&c| c as usize).sum::<usize>() < self.threshold2 {
        //                 pref -= 1;
        //             } else {
        //                 break;
        //             }
        //         }
        //     }
        // }

        // if pref == 1 {
        // right-max kmer?
        // if (kmer[self.k - 1] as usize) < self.threshold {
        //     pref = 2;
        // }
        // for &tmer in &kmer[..self.k - 1] {
        //     if (tmer as usize) >= self.threshold {
        //         pref = 2;
        //     }
        // }
        // }
        // let v = kmer[0] as isize - *kmer[1..].iter().max().unwrap() as isize;
        // } else {
        //     pref = 1;
        // }
        //
        // Classic max-in-middle open variant.
        //
        // if self.open && pref as usize >= 100 {
        //     let m = self.k / 2;
        //     for (i, &tmer) in kmer.iter().enumerate() {
        //         if ((tmer as usize) >= self.threshold) ^ (i == m) {
        //             pref = 200;
        //         }
        //     }
        // }
        (
            pref,
            0,
            // -v,
            // (255 - kmer[0]) as isize,
            // kmer[0] as isize - m as isize,
            Order::key(&self.o, kmer),
        )
    }
}

impl SamplingScheme for ThresholdMinimizer {
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let mut good = false;

        let h = self.hash_kmer(&lmer[0..self.k]);
        if h.0 < 100 {
            self.count3.set(self.count3.get() + 1);
        }

        // Do not allow kmers that are followed by only 0s, since they are likely to overlap a 10000.
        let mut skip = 0;
        if self.open {}

        let ((pref, _, _), mut p, _) = lmer
            .windows(self.k)
            .enumerate()
            .map(|(i, kmer)| {
                let h = self.hash_kmer(kmer);
                if h.0 < 100 {
                    self.count.set(self.count.get() + 1);
                    good = true;
                }
                (h, i, kmer)
            })
            .min()
            .unwrap();

        // if self.open && pref == 100 && self.k <= self.w {
        //     // Check if maybe there is a large char in the last k-1 positions.
        //     // Then take the kmer ending just before it.
        //     for (i, &c) in lmer[self.w..].iter().enumerate() {
        //         if c >= self.threshold as u8 {
        //             p = i + self.w - self.k;
        //         }
        //     }
        // }

        if good {
            self.count2.set(self.count2.get() + 1);
        }
        p
    }
}

impl Drop for ThresholdMinimizer {
    fn drop(&mut self) {
        eprintln!("k {} => count : {}", self.k, self.count.get());
        eprintln!("k {} => count2: {}", self.k, self.count2.get());
        eprintln!("k {} => count3: {}", self.k, self.count3.get());
    }
}
