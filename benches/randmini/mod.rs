pub mod baseline;
pub mod daniel;
pub mod queue;
pub mod rescan;

pub trait Minimizer {
    fn minimizers(&self, text: &[u8]) -> Vec<usize>;
}

pub trait SlidingMin<V> {
    /// Initialize a new datastructure with window size `w`.
    fn new(w: usize) -> Self;
    /// Push a new value, starting at position 0.
    /// Return the pos and value of the minimum of the last w elements.
    fn push(&mut self, val: V) -> Elem<V>;
}

#[derive(Clone, Copy)]
pub struct Elem<Val> {
    pub pos: usize,
    pub val: Val,
}

trait Max {
    const MAX: Self;
}

impl Max for usize {
    const MAX: usize = usize::MAX;
}
