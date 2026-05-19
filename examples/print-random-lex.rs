use itertools::Itertools;
use minimizers::order::{Order, RandomLexO};

fn main() {
    let args = std::env::args().collect_vec();
    let seed = args[1].parse::<usize>().unwrap();
    let k = args[2].parse::<usize>().unwrap();
    let o = RandomLexO(seed as u128);

    let sigma = 2;
    let mut strings = vec![];
    for k in 1..=k {
        let kmer_it = std::iter::repeat(0..sigma)
            .take(k)
            .multi_cartesian_product();
        'kmer: for kmer in kmer_it {
            // filter kmers containing a smaller suffix
            let h = o.key(&kmer);
            for i in 1..k {
                let suffix = &kmer[i..];
                if o.key(suffix) < h {
                    continue 'kmer;
                }
            }

            strings.push(kmer);
        }
    }
    strings.sort_by_cached_key(|s| o.key(&s));
    for x in strings {
        for c in x {
            eprint!("{}", c);
        }
        eprintln!();
    }
}
