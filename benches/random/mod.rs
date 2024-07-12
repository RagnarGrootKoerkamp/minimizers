mod baseline;
pub use baseline::*;

pub trait Minimizer {
    fn minimizers(&self, text: &[u8]) -> Vec<usize>;
}
