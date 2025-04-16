//! Take a .fa file, compute minimizers, and highlight those in red.

use std::{ops::Range, thread::sleep, time::Duration};

use clap::{Parser, Subcommand};
use colored::Colorize;
use minimizers::{order::Lex, schemes::Minimizer, SamplingScheme};
use needletail::Sequence;

#[derive(Parser)]
enum Commands {
    /// Adds files to myapp
    Mini { input: String, k: usize, w: usize },
    Print {
        input: String,
        #[clap(long)]
        highlight: bool,
        #[clap(long)]
        virus: bool,
        #[clap(long, default_value_t = 0)]
        skip: usize,
    },
}

fn main() {
    let args = Commands::parse();

    match args {
        Commands::Mini { .. } => mini(args),
        Commands::Print {
            input,
            highlight,
            virus,
            skip,
        } => {
            let mut reader = needletail::parse_fastx_file(input).unwrap();
            let record = reader.next().unwrap().unwrap();
            let text = record.seq();
            let text = &text[skip..];

            if !highlight {
                print!("{}\n", as_str(&text).black().on_white());
            } else {
                for &c in text.iter() {
                    if virus {
                        match c {
                            b'A' => print!("{}", as_str(&[c]).truecolor(0, 255, 0).on_black()),
                            b'C' => print!("{}", as_str(&[c]).truecolor(0, 208, 0).on_black()),
                            b'G' => print!("{}", as_str(&[c]).truecolor(0, 155, 0).on_black()),
                            b'T' => print!("{}", as_str(&[c]).white().on_black()),
                            _ => panic!(),
                        }
                    } else {
                        match c {
                            b'A' => print!(
                                "{}",
                                as_str(&[c])
                                    .truecolor(255, 0, 0)
                                    .on_truecolor(255, 255, 255)
                            ),
                            b'C' => print!(
                                "{}",
                                as_str(&[c])
                                    .truecolor(0, 128, 0)
                                    .on_truecolor(255, 255, 255)
                            ),
                            b'G' => print!(
                                "{}",
                                as_str(&[c])
                                    .truecolor(0, 0, 255)
                                    .on_truecolor(255, 255, 255)
                            ),
                            b'T' => print!(
                                "{}",
                                as_str(&[c]).truecolor(0, 0, 0).on_truecolor(255, 255, 255)
                            ),
                            _ => panic!(),
                        }
                    }
                }
                println!();
            }
        }
    }
}

fn mini(args: Commands) {
    let Commands::Mini { input, k, w } = args else {
        panic!();
    };

    let mut reader = needletail::parse_fastx_file(input).unwrap();
    let record = reader.next().unwrap().unwrap();
    let text = &record.seq();

    highlight(text, w, k);
}

fn highlight(text: &[u8], w: usize, k: usize) {
    let m = Minimizer::new(k, w, Lex);
    let poss = m.stream(text);
    let ranges = poss
        .into_iter()
        .map(|start| start..start + k)
        .collect::<Vec<_>>();
    // Join overlapping/touching ranges.
    let ranges = ranges
        .into_iter()
        .fold(Vec::new(), |mut acc: Vec<Range<usize>>, range| {
            if let Some(last) = acc.last_mut() {
                if last.end >= range.start {
                    last.end = last.end.max(range.end);
                } else {
                    acc.push(range);
                }
            } else {
                acc.push(range);
            }
            acc
        });

    // Print the sequence with minimizers highlighted.
    let mut start = 0;
    for range in ranges {
        if start < range.start {
            print!(
                "{}",
                as_str(&text[start..range.start])
                    .black()
                    .on_truecolor(255, 255, 255)
            );
        }
        print!("{}", as_str(&text[range.clone()]).black().on_yellow());
        start = range.end;
    }
    print!(
        "{}",
        as_str(&text[start..]).black().on_truecolor(255, 255, 255)
    );
    print!(
        "{}",
        format!("\nk = {k}   w = {w}\n")
            .black()
            .on_truecolor(255, 255, 255)
    );
}

fn as_str(s: &[u8]) -> String {
    String::from_utf8_lossy(s).to_string()
}
