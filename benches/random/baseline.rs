use super::*;

pub struct V0NaiveLex {
    pub w: usize,
    pub k: usize,
}

impl Minimizer for V0NaiveLex {
    fn minimizers(&self, text: &[u8]) -> Vec<usize> {
        // Iterate over the windows of size l=w+k-1.
        text.windows(self.w + self.k - 1)
            .enumerate()
            // For each window, starting at pos j, find the lexicographically smallest k-mer.
            .map(|(j, window)| {
                j + window
                    .windows(self.k)
                    .enumerate()
                    .min_by_key(|(_idx, kmer)| *kmer)
                    .unwrap()
                    .0
            })
            .collect()
    }
}

pub struct V1NaiveFx {
    pub w: usize,
    pub k: usize,
}

impl Minimizer for V1NaiveFx {
    fn minimizers(&self, text: &[u8]) -> Vec<usize> {
        // Iterate over the windows of size l=w+k-1.
        text.windows(self.w + self.k - 1)
            .enumerate()
            // For each window, starting at pos j, find the lexicographically smallest k-mer.
            .map(|(j, window)| {
                j + window
                    .windows(self.k)
                    .enumerate()
                    .min_by_key(|(_idx, kmer)| fxhash::hash64(kmer))
                    .unwrap()
                    .0
            })
            .collect()
    }
}

pub struct V2NaiveWy {
    pub w: usize,
    pub k: usize,
}

impl Minimizer for V2NaiveWy {
    fn minimizers(&self, text: &[u8]) -> Vec<usize> {
        // Iterate over the windows of size l=w+k-1.
        text.windows(self.w + self.k - 1)
            .enumerate()
            // For each window, starting at pos j, find the lexicographically smallest k-mer.
            .map(|(j, window)| {
                j + window
                    .windows(self.k)
                    .enumerate()
                    .min_by_key(|(_idx, kmer)| wyhash::wyhash(kmer, 0))
                    .unwrap()
                    .0
            })
            .collect()
    }
}
