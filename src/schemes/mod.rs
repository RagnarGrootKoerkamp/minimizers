#![allow(refining_impl_trait)]

mod anchors;
mod bruteforce;
mod decycling;
mod explicit;
mod fracmin;
mod minimizer;
mod modmini;
mod open_closed;
mod rot_mini;
mod threshold;

pub use anchors::*;
pub use bruteforce::*;
pub use decycling::*;
pub use explicit::*;
pub use fracmin::*;
pub use minimizer::*;
pub use modmini::*;
pub use open_closed::*;
pub use rot_mini::*;
pub use threshold::*;

use super::{order::*, MonotoneQueue, Params, SamplingScheme};
use serde::{Deserialize, Serialize};

/// Classic random minimizers.
pub type RandomMinimizer = minimizer::Minimizer<RandomO>;

// GENERIC MINIMIZER SCHEME.

/// Generic minimizer scheme.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct M<O: ToOrder>(pub O);
/// Wrapper that does tiebreaking by kmer hash as last step.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct RM<O: ToOrder>(pub O);

#[typetag::serialize]
impl<O: ToOrder> Params for M<O> {
    fn build(&self, w: usize, k: usize, sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(Minimizer::new(k, w, self.0.to_order(w, k, sigma)))
    }
}
impl<O: ToOrder> M<O> {
    fn build_from_order(o: &O, w: usize, k: usize, sigma: usize) -> Minimizer<O::O> {
        Minimizer::new(k, w, o.to_order(w, k, sigma))
    }
}
#[typetag::serialize]
impl<O: ToOrder> Params for RM<O> {
    fn build(&self, w: usize, k: usize, sigma: usize) -> Box<dyn SamplingScheme> {
        M((self.0.clone(), RandomO)).build(w, k, sigma)
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
