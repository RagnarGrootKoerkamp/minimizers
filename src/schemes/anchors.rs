use super::*;

#[typetag::serialize]
impl Params for BdAnchor {
    fn build(&self, w: usize, k: usize, _sigma: usize) -> Box<dyn SamplingScheme> {
        let r = self.r.max(k - 1);
        // assert!(self.r >= k - 1);
        Box::new(BdAnchorS::new(w + k - 1, r))
    }
}

#[typetag::serialize]
impl<O: ToOrder> Params for SusAnchor<O> {
    fn build(&self, w: usize, k: usize, sigma: usize) -> Box<dyn SamplingScheme> {
        Box::new(SusAnchorS::new(w, k, self.0.to_order(w, k, sigma)))
    }
}

pub struct BdAnchorS {
    w: usize,
    r: usize,
}

impl BdAnchorS {
    pub fn new(w: usize, r: usize) -> Self {
        Self { w, r }
    }
}

impl SamplingScheme for BdAnchorS {
    fn l(&self) -> usize {
        self.w
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let w = self.w;
        debug_assert_eq!(lmer.len(), w);
        let mut best = 0;
        for i in 1..w.saturating_sub(self.r) {
            for j in 0..w {
                if lmer[(i + j) % w] != lmer[(best + j) % w] {
                    if lmer[(i + j) % w] < lmer[(best + j) % w] {
                        best = i;
                    }
                    break;
                }
            }
        }
        best
    }
}

/// NOTE: O should be Lex or AntiLex order. Random order will not be good.
pub struct SusAnchorS<O: Order> {
    w: usize,
    k: usize,
    o: O,
}

impl<O: Order> SusAnchorS<O> {
    pub fn new(w: usize, k: usize, o: O) -> Self {
        Self { w, k, o }
    }
}

impl<O: Order> SamplingScheme for SusAnchorS<O> {
    fn l(&self) -> usize {
        self.w + self.k - 1
    }

    fn sample(&self, lmer: &[u8]) -> usize {
        let mut best = (self.o.key(lmer), 0);
        for i in 1..self.w {
            best = best.min((self.o.key(&lmer[i..]), i));
        }
        best.1
    }

    fn stream(&self, text: &[u8]) -> Vec<usize> {
        self.stream_queue(text)
    }
}

impl<O: Order> SusAnchorS<O> {
    pub fn stream_events(&self, text: &[u8]) -> Vec<usize> {
        let w = self.w;
        let k = self.k;
        let l = self.l();
        let mut active = vec![];
        let events_len = 2 * self.l();
        let mut events = vec![vec![]; events_len];

        let mut poss = vec![];

        for i in 0..text.len() - (k - 1) {
            let end = i + k;
            let mut evs = std::mem::take(&mut events[(end - 1) % events_len]);
            // Process the new empty suffix.
            active.push(i);
            evs.push(i);
            // eprintln!("push {i}: {:?}", &text[i..(i + l).min(text.len())]);
            // eprintln!("active: {active:?}");
            // eprintln!("events: {evs:?}");
            // Process events
            for &pos in &evs {
                if let Some(mut idx) = active.iter().position(|x| *x == pos) {
                    while idx > 0 {
                        let prev_active = active[idx - 1];
                        // eprintln!("cmp prev {prev_active}.. with new {pos}..");
                        let new_suffix = &text[pos..end];
                        // eprintln!("new ..: {new_suffix:?}");
                        let old_suffix = &text[prev_active..end];
                        // eprintln!("prev..: {old_suffix:?}");
                        if self.o.key(new_suffix) < self.o.key(old_suffix) {
                            // eprintln!("new is smaller; drop prev");
                            active.remove(idx - 1);
                            idx -= 1;
                        } else {
                            // Check if new suffix becomes smaller eventually
                            let l = lcp(
                                &text[pos..(pos + l).min(text.len())],
                                &text[prev_active..(prev_active + l).min(text.len())],
                            );
                            // eprintln!("new is not larger. lcp={l}");
                            if l < self.l()
                                && pos + l < text.len()
                                && self.o.key(&text[pos..=pos + l])
                                    < self.o.key(&text[prev_active..=prev_active + l])
                            {
                                // eprintln!("new will be smaller from i={}", pos + l);
                                events[(pos + l) % events_len].push(pos);
                            }
                            break;
                        }
                    }
                }
            }
            evs.clear();
            events[(end - 1) % events_len] = evs;

            // eprintln!("{i}: {active:?}");
            assert!(i.saturating_sub(w - 1) <= active[0] && active[0] <= i);
            if i >= w - 1 {
                poss.push(active[0]);
            }
            if active[0] + (w - 1) == i {
                active.remove(0);
            }
        }
        poss
    }
    pub fn stream_queue(&self, text: &[u8]) -> Vec<usize> {
        let w = self.w;
        let k = self.k;
        let l = self.l();
        // (start-pos of anchor, index of text character where it becomes active)
        let mut active: Vec<(usize, usize)> = vec![];

        let mut poss = vec![];

        for i in 0..text.len() - (k - 1) {
            let end = i + k;
            // Process the new empty suffix.
            // eprintln!("push {i}: {:?}", &text[i..(i + l).min(text.len())]);
            // eprintln!("active: {active:?}");
            {
                // Find out where this suffix becomes active
                let pos = i;
                let mut active_at = end - 1;
                while let Some(&(prev, prev_activates_at)) = active.last() {
                    // eprintln!("cmp prev {prev} ({prev_activates_at}).. with new {pos}..");
                    // let new_suffix = &text[pos..end];
                    // eprintln!("new ..: {new_suffix:?}");
                    // let old_suffix = &text[prev..end];
                    // eprintln!("prev..: {old_suffix:?}");
                    // Check if new suffix becomes smaller eventually
                    let l = lcp(
                        &text[pos..(pos + l).min(text.len())],
                        &text[prev..(prev + l).min(text.len())],
                    );
                    // eprintln!("lcp={l}");
                    if pos + l < text.len()
                        && self.o.key(&text[pos..=pos + l]) < self.o.key(&text[prev..=prev + l])
                    {
                        let pos_activates_at = (pos + l).min(prev + self.l());
                        // assert!(pos_activates_at != prev_activates_at);
                        if pos_activates_at >= end && pos_activates_at > prev_activates_at {
                            // eprintln!(
                            // "New will be smaller later: active from {pos_activates_at} > {prev_activates_at}"
                            // );
                            active_at = pos_activates_at;
                            break;
                        } else {
                            // eprintln!(
                            // "New will be smaller before: active from {pos_activates_at} < {prev_activates_at}. Pop prev."
                            // );
                            active.pop();
                            if active.is_empty() {
                                // active_at = end - 1;
                                break;
                            }
                            continue;
                        }
                    } else {
                        // eprintln!("New is larger: active when previous element drops");
                        active_at = prev + self.l();
                        break;
                    }
                    // unreachable
                }
                // eprintln!("Pos {pos} activates at char {active_at}");
                assert!(pos <= active_at);
                active.push((pos, active_at));
            }

            // eprintln!("{i}: {active:?}");
            assert!(i.saturating_sub(w - 1) <= active[0].0 && active[0].0 <= i);
            if i >= w - 1 {
                poss.push(active[0].0);
            }
            if active[0].0 + (w - 1) == i {
                active.remove(0);
            } else if active.len() > 1 && active[1].1 <= (i + k) {
                active.remove(0);
            }
        }
        poss
    }
}

// TODO: Word-based or SIMD-based
fn lcp(a: &[u8], b: &[u8]) -> usize {
    let l = a.len().min(b.len());
    for i in 0..l {
        if a[i] != b[i] {
            return i;
        }
    }
    return l;
}
