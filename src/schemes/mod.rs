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

#[derive(Debug, Clone, Serialize)]
pub struct OpenClosed<O: ToOrder> {
    pub r: usize,
    /// Do open syncmers first?
    pub open: bool,
    /// Do closed syncmers second?
    pub closed: bool,
    /// When true, open syncmers must have this offset instead of in the middle by default.
    pub offset: Option<usize>,
    /// When true, any position offset%w makes a kmer an open syncmer.
    pub modulo: bool,
    pub open_by_tmer: bool,
    pub closed_by_tmer: bool,
    pub other_by_tmer: bool,
    pub anti_tmer: bool,
    pub o: O,
    /// r = max(r, k-w) if true.
    pub miniception_r: bool,
}

/// Asymptotic (in k) optimal minimizers:
/// - Assume k=x*w; sum all i mod w positions.
/// - Take minimizer if first coordinate is max or at most sigma-1 away from max.
/// - If multiple, take random.
///
/// Assumes alphabet size sigma=4.
/// Also requires the input to be encoded as values 0..sigma, not as arbitrary ascii codepoints.
///
/// On the overshoot:
/// Instead of the original sigma, we can use sigma-1, since the max difference
/// between two characters is sigma-1.  In fact, I conjecture that
/// floor((sigma+1)/2) overshoot is sufficient. I have a sketch of a proof for
/// this.
#[derive(Clone, Serialize, Debug)]
pub struct RotMinimizer;

pub struct RotMinimizerO {
    w: usize,
    sigma: usize,
}

// ANCHORS.

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BdAnchor {
    pub r: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SusAnchor<O: ToOrder>(pub O);

// MINIMIZER WRAPPERS.

#[derive(Debug, Serialize)]
pub struct ModP {
    pub r: usize,
    pub lr: bool,
    pub t: usize,
    pub params: Box<dyn Params>,
}
