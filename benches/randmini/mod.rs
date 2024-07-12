pub mod baseline;
pub mod daniel;
pub mod sliding_min;

pub trait Minimizer {
    fn minimizers(&self, text: &[u8]) -> Vec<usize>;
}
