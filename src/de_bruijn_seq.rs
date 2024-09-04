/// A function to generate a De Bruijn sequence.
/// Transcribed from wikipedia; I don't really understand the code.
/// https://en.wikipedia.org/wiki/De_Bruijn_sequence
///
/// sigma: alphabet size
/// order: order of the De Bruijn sequence. All sigma^order strings occur once.
/// The returned string has length sigma^order + order - 1.
pub fn de_bruijn_sequence(sigma: usize, order: usize) -> Vec<u8> {
    let mut a = vec![0; order * sigma];
    let mut s = Vec::with_capacity(sigma.pow(order as u32) + order);

    fn recursion(sigma: usize, order: usize, t: usize, p: usize, a: &mut Vec<u8>, s: &mut Vec<u8>) {
        if t > order {
            if order % p == 0 {
                s.extend(&a[1..p + 1]);
            }
        } else {
            a[t] = a[t - p];
            recursion(sigma, order, t + 1, p, a, s);
            for c in a[t - p] + 1..sigma as u8 {
                a[t] = c;
                recursion(sigma, order, t + 1, t, a, s);
            }
        }
    }

    recursion(sigma, order, 1, 1, &mut a, &mut s);
    s.extend_from_within(0..order - 1);
    s
}

/// Returns a string that measures the exact density.
/// The returned string has length `sigma^order + (k-1)`, so that it contains `sigma^order` kmers.
///
/// Following Marcais' 2017 and 2018:
/// For forward schemes: De Bruijn seq of order w+k.
/// For local schemes: De Bruijn seq of order 2w+k-2.
pub fn exact_density_string(k: usize, w: usize, sigma: usize, forward: bool) -> Vec<u8> {
    // FIXME +1?
    let order = if forward { w + k } else { 2 * w + k - 2 } + 1;
    let mut s = de_bruijn_sequence(sigma, order);
    s.truncate(s.len() - (order - 1) + (k - 1));
    s
}

/// Returns a string of sufficient length to measure the exact density.
/// Also returns the number of lmers that should be considered.
pub fn cyclic_exact_density_string(
    k: usize,
    w: usize,
    sigma: usize,
    forward: bool,
) -> (Vec<u8>, usize) {
    let order = if forward { w + k } else { 2 * w + k - 2 };
    let s = de_bruijn_sequence(sigma, order);
    (s, sigma.pow(order as u32))
}

mod test {
    #[test]
    fn test_de_bruin_sequence() {
        let s = super::de_bruijn_sequence(2, 3);
        assert_eq!(s, [0, 0, 0, 1, 0, 1, 1, 1, 0, 0]);
        let s = super::de_bruijn_sequence(3, 2);
        assert_eq!(s, [0, 0, 1, 0, 2, 1, 1, 2, 2, 0]);
    }
}
