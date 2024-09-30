#![allow(refining_impl_trait)]

mod anchors;
mod bruteforce;
mod decycling;
mod explicit;
mod fracmin;
mod miniception;
mod minimizer;
mod modmini;
mod rot_mini;
mod syncmer;
mod threshold;

pub use anchors::*;
pub use bruteforce::*;
pub use decycling::*;
pub use explicit::*;
pub use fracmin::*;
pub use miniception::*;
pub use minimizer::*;
pub use modmini::*;
pub use rot_mini::*;
pub use syncmer::*;
pub use threshold::*;

use super::{order::*, MonotoneQueue, Params, SamplingScheme};
use serde::{Deserialize, Serialize};

/// Classic random minimizers.
pub type RandomMinimizer = minimizer::Minimizer<RandomO>;

// GENERIC MINIMIZER SCHEME.

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct M<O: ToOrder>(pub O);

#[typetag::serialize]
impl<O: ToOrder + Serialize + 'static> Params for M<O> {
    fn build(&self, w: usize, k: usize, sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(Minimizer::new(k, w, self.0.to_order(w, k, sigma)))
    }
}

// ORDERS.

/// (Double) decycling minimizers based of the Mykkeltveit embedding.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Decycling {
    pub double: bool,
}

/// Sample 1/f of the kmers.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FracMin {
    pub f: usize,
}

// MINIMIZER WRAPPERS.

#[derive(Debug, Serialize)]
pub struct ModP {
    pub r: usize,
    pub lr: bool,
    pub t: usize,
    pub params: Box<dyn Params>,
}
