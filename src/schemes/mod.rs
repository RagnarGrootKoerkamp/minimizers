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
pub type RandomMinimizer = minimizer::Minimizer<RandomOrder>;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RandomM;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AntiLexM;

#[typetag::serde]
impl Params for RandomM {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(Minimizer::new(k, w, RandomOrder))
    }
}

#[typetag::serde]
impl Params for AntiLexM {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(Minimizer::new(k, w, AntiLex))
    }
}
