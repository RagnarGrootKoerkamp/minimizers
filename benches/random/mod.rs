pub mod baseline;
pub mod daniel;

pub trait Minimizer {
    fn minimizers(&self, text: &[u8]) -> Vec<usize>;
}
