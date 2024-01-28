eval:
    cargo run -r -- -n 300000 -s 4 eval -o data/results_4.json
    cargo run -r -- -n 300000 -s 256 eval -o data/results_256.json

run:
    cargo run -r -- -n 10000 -s 256 eval

flame:
    cargo flamegraph --open --skip-after minimizers::MinimizerType::density -- -n 10000 -s 256 eval
