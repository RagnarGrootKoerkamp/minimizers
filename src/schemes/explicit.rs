use super::*;

#[derive(Debug, Clone)]
pub struct ExplicitLocalScheme {
    pub k: usize,
    pub w: usize,
    pub sigma: usize,
    /// The index in [w] to choose for each of sigma^(k+w-1) possible l-mers.
    pub map: Vec<u8>,
}

impl SamplingScheme for ExplicitLocalScheme {
    fn w(&self) -> usize {
        self.w
    }
    fn k(&self) -> usize {
        self.k
    }
    fn l(&self) -> usize {
        self.k + self.w - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        debug_assert_eq!(lmer.len(), self.l());
        debug_assert_eq!(self.map.len(), self.sigma.pow(self.l() as u32));
        let mut v = 0;
        for c in lmer {
            assert!(*c < self.sigma as u8);
            v = self.sigma * v + *c as usize;
        }
        self.map[v] as usize
    }
}

impl ExplicitLocalScheme {
    /// Check if the local scheme corresponds to a directed order on kmers.
    /// If so, return one such order.
    pub fn to_directed_order(&self) -> Option<ExplicitDirectedOrder> {
        let l = self.k + self.w - 1;
        let num_lmers = self.sigma.pow(l as u32);
        assert_eq!(self.map.len(), num_lmers);
        let num_kmers = self.sigma.pow(self.k as u32);
        // smaller[i] is a list of kmers j with i < j, indicating that i comes before j in the toposort.
        let mut smaller = vec![vec![]; num_kmers];
        // direction in which kmer i should be compared in case of equality.
        let mut directions = vec![None; num_kmers];

        for lmer in 0..num_lmers {
            let mut chars = vec![];
            {
                let mut lmer = lmer;
                for _ in 0..l {
                    chars.push((lmer % self.sigma) as u8);
                    lmer /= self.sigma;
                }
                chars.reverse();
            }

            let kmers = {
                let mut lmer = lmer;
                let mut kmers = (0..self.w)
                    .map(|_| {
                        let kmer = lmer % num_kmers;
                        lmer /= self.sigma;
                        kmer
                    })
                    .collect::<Vec<_>>();
                kmers.reverse();
                kmers
            };

            let chosen_idx = self.map[lmer] as usize;
            let chosen_kmer = kmers[chosen_idx];
            // eprintln!("lmer {lmer:>3}: {:?} => {chosen_idx}", chars);
            for i in 0..self.w {
                if i == chosen_idx {
                    continue;
                }
                let kmer = kmers[i];
                if kmer != chosen_kmer {
                    smaller[chosen_kmer].push(kmer);
                    // eprintln!(
                    //     "{:?} < {:?}",
                    //     &chars[chosen_idx..chosen_idx + self.k],
                    //     &chars[i..i + self.k]
                    // );
                } else {
                    let wanted_direction = if i < chosen_idx {
                        Some(Direction::Rightmost)
                    } else {
                        Some(Direction::Leftmost)
                    };
                    if directions[chosen_kmer].is_none() {
                        directions[chosen_kmer] = wanted_direction;
                    } else {
                        if directions[chosen_kmer] != wanted_direction {
                            eprintln!(
                                "Inconsistent direction for kmer {:?}",
                                &chars[i..i + self.k]
                            );
                            // return None;
                        }
                    }
                }
            }
        }

        for x in &mut smaller {
            x.sort();
            x.dedup();
        }

        // Toposort the 'smaller' graph.
        let mut order = vec![];
        {
            let mut visited = vec![0; num_kmers];
            fn dfs(
                i: usize,
                smaller: &Vec<Vec<usize>>,
                visited: &mut Vec<u8>,
                order: &mut Vec<usize>,
            ) -> Option<()> {
                if visited[i] == 1 {
                    eprintln!("Cycle detected");
                    return None;
                }
                if visited[i] == 2 {
                    return Some(());
                }
                visited[i] = 1;
                for &j in &smaller[i] {
                    dfs(j, smaller, visited, order)?;
                }
                visited[i] = 2;
                order.push(i);
                Some(())
            }
            for i in 0..num_kmers {
                dfs(i, &smaller, &mut visited, &mut order)?;
            }
            order.reverse();
        }

        let mut idx = vec![(0, Direction::Leftmost); num_kmers];
        for (i, &x) in order.iter().enumerate() {
            idx[x].0 = i;
            idx[x].1 = directions[x].unwrap_or(Direction::Leftmost);
        }

        Some(ExplicitDirectedOrder {
            k: self.k,
            sigma: self.sigma,
            idx,
        })
    }
}
