use itertools::Itertools;
use std::marker::PhantomData;

pub mod minimizer;
pub mod naive;
pub mod sliding_min;

use minimizer::*;
use naive::*;
use sliding_min::*;

trait Max {
    const MAX: Self;
}

impl Max for usize {
    const MAX: usize = usize::MAX;
}
impl Max for u64 {
    const MAX: u64 = u64::MAX;
}
