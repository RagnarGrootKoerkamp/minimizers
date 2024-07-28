use itertools::Itertools;
use std::marker::PhantomData;

pub mod counting;
pub mod hash;
pub mod jumping;
pub mod minimizer;
pub mod naive;
pub mod queue;
pub mod rescan;
pub mod rescan2;
pub mod ringbuf;
pub mod sliding_min;
pub mod split;

pub use counting::*;
pub use hash::*;
pub use jumping::*;
pub use minimizer::*;
pub use naive::*;
pub use queue::*;
pub use rescan::*;
pub use rescan2::*;
pub use ringbuf::*;
pub use sliding_min::*;
pub use split::*;

trait Max {
    const MAX: Self;
}

impl Max for usize {
    const MAX: usize = usize::MAX;
}
impl Max for u64 {
    const MAX: u64 = u64::MAX;
}
