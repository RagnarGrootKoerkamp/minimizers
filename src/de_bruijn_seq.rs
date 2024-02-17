/// A function to generate a De Bruijn sequence.
/// Transcribed from wikipedia; I don't really understand the code.
/// https://en.wikipedia.org/wiki/De_Bruijn_sequence
///
/// sigma: alphabet size
/// order: order of the De Bruijn sequence. All sigma^order strings occur once.
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

mod test {
    #[test]
    fn test_de_bruin_sequence() {
        let s = super::de_bruijn_sequence(2, 3);
        assert_eq!(s, [0, 0, 0, 1, 0, 1, 1, 1]);
        let s = super::de_bruijn_sequence(3, 2);
        assert_eq!(s, [0, 0, 1, 0, 2, 1, 1, 2, 2]);
    }
}
