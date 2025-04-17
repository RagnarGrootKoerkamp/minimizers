use std::{fs::File, io::Read};

use serde::Serialize;

use crate::order::{Order, ToOrder};

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct GreedyP;

pub struct Greedy {
    base: usize,
    factor: usize,
    idx: Vec<u16>,
}

impl ToOrder for GreedyP {
    type O = Greedy;

    fn to_order(&self, w: usize, k: usize, _sigma: usize) -> Self::O {
        let paths = [
            format!("greedy-mini-data/w{w}_k{k}.bin"),
            format!("../greedy-mini-data/w{w}_k{k}.bin"),
        ];
        let mut data = vec![];
        'data: {
            for p in paths {
                let f = File::open(p);
                if f.is_err() {
                    continue;
                }
                f.unwrap().read_to_end(&mut data).unwrap();
                break 'data;
            }
            panic!("greedy-mini-data/w{w}_k{k}.bin not found!\nRun from root directory or one of its direct subdirectories.");
        }
        let len = usize::from_ne_bytes(data[0..8].try_into().unwrap());
        assert_eq!(data.len(), 8 * (len + 1));
        let mut idx = vec![0; len];
        let base = 4usize;
        for i in 0..len {
            let start = 8 * (i + 1);
            let end = start + 8;
            idx[i] = u64::from_ne_bytes(data[start..end].try_into().unwrap()) as u16;
        }
        let factor = base.pow(k as _);
        Greedy { base, factor, idx }
    }
}

impl Order for Greedy {
    type T = usize;

    fn key(&self, kmer: &[u8]) -> Self::T {
        // Order on sigma=2 lower-bit alphabet, via greedy scheme.
        let mut x = 0;
        // Order on remaining bits, lex.
        let mut y = 0;
        for &c in kmer {
            y = y * self.base + ((c >> 1) as usize & (self.base - 1));
            x = (x << 1) | (c & 1) as usize;
        }
        self.factor * self.idx[x] as usize + y
    }
}
