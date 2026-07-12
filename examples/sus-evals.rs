//! Speed and density evals for SUS-anchors on random and human data.

use minimizers::{
    generate_true_random_string,
    order::{AntiLex, Lex, RandomO, ToOrder},
    schemes::{greedy::GreedyP, BdAnchorS, Decycling, Minimizer, SusAnchorS},
    SamplingScheme,
};

fn main() {
    let len;
    let sigma = 4;

    let hg = {
        let mut reader = needletail::parse_fastx_file("human-genome.fa").unwrap();
        let mut text = vec![];
        while let Some(record) = reader.next() {
            let record = record.unwrap();
            let id = str::from_utf8(record.id()).unwrap();
            if id.contains("chromosome 19") {
                eprintln!("{:?}", id);
                text.extend_from_slice(&record.seq());
            }
        }
        len = text.len();
        text
    };
    let text = generate_true_random_string(len, sigma);
    eprintln!("len {len}");

    for w in [31, 63] {
        eprintln!("w = {w}");
        let ss: Vec<(&str, Box<dyn SamplingScheme>)> = vec![
            // ("sus-lex", Box::new(SusAnchorS::new(w, 1, Lex))),
            // ("sus-antilex", Box::new(SusAnchorS::new(w, 1, AntiLex))),
            // ("bd-2", Box::new(BdAnchorS::new(w - 1, 2))),
            // ("bd-3", Box::new(BdAnchorS::new(w - 2, 3))),
            // ("bd-4", Box::new(BdAnchorS::new(w - 3, 4))),
            // // ("sus-abb", Box::new(SusAnchorS::new(w, 1, ABB))),
            // // ("sus-abb+", Box::new(SusAnchorS::new(w, 1, (ABB, Lex)))),
            // ("min-lex-2", Box::new(Minimizer::new(2, w - 1, Lex))),
            // ("min-lex-3", Box::new(Minimizer::new(3, w - 2, Lex))),
            // ("min-lex-4", Box::new(Minimizer::new(4, w - 3, Lex))),
            // ("min-lex-5", Box::new(Minimizer::new(5, w - 4, Lex))),
            // ("min-antilex-2", Box::new(Minimizer::new(2, w - 1, AntiLex))),
            // ("min-antilex-3", Box::new(Minimizer::new(3, w - 2, AntiLex))),
            // ("min-antilex-4", Box::new(Minimizer::new(4, w - 3, AntiLex))),
            // ("min-antilex-5", Box::new(Minimizer::new(5, w - 4, AntiLex))),
            // ("min-random-2", Box::new(Minimizer::new(2, w - 1, RandomO))),
            // ("min-random-3", Box::new(Minimizer::new(3, w - 2, RandomO))),
            // ("min-random-4", Box::new(Minimizer::new(4, w - 3, RandomO))),
            // ("min-random-5", Box::new(Minimizer::new(5, w - 4, RandomO))),
            // (
            //     "double-decycling-3",
            //     Box::new(Minimizer::new(
            //         3,
            //         w - 2,
            //         (
            //             Decycling { double: true }.to_order(w - 2, 3, sigma),
            //             RandomO,
            //         ),
            //     )),
            // ),
            // (
            //     "double-decycling-4",
            //     Box::new(Minimizer::new(
            //         4,
            //         w - 3,
            //         (
            //             Decycling { double: true }.to_order(w - 3, 4, sigma),
            //             RandomO,
            //         ),
            //     )),
            // ),
            // (
            //     "double-decycling-5",
            //     Box::new(Minimizer::new(
            //         5,
            //         w - 4,
            //         (
            //             Decycling { double: true }.to_order(w - 4, 5, sigma),
            //             RandomO,
            //         ),
            //     )),
            // ),
            (
                "greedy-7",
                Box::new(Minimizer::new(7, w - 6, GreedyP.to_order(w - 6, 7, sigma))),
            ),
            (
                "greedy-6",
                Box::new(Minimizer::new(6, w - 5, GreedyP.to_order(w - 5, 6, sigma))),
            ),
            (
                "greedy-5",
                Box::new(Minimizer::new(5, w - 4, GreedyP.to_order(w - 4, 5, sigma))),
            ),
            (
                "greedy-4",
                Box::new(Minimizer::new(4, w - 3, GreedyP.to_order(w - 3, 4, sigma))),
            ),
            (
                "greedy-3",
                Box::new(Minimizer::new(3, w - 2, GreedyP.to_order(w - 2, 3, sigma))),
            ),
        ];
        for (name, s) in ss {
            let [(d1, t1), (d2, _t2)] = [&text, &hg].map(|text| {
                let start = std::time::Instant::now();
                let sampled = s.stream(&text);
                let end = std::time::Instant::now();
                let c = 1 + sampled.array_windows().filter(|[a, b]| a != b).count();
                let density = c as f32 / len as f32 * (w + 1) as f32;

                let t = len as f32 / (end - start).as_secs_f32() / 1e6;
                (density, t)
            });
            let lb = 2. / (w + 1) as f32;
            eprintln!(
                "{name:<20}: lb {lb:.6} density: {d1:.4} {d2:.4}, throughput: {t1:.1?} Mbp/s",
            );
        }
    }
}
