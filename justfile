eval:
    # cargo run -r -- -n 1000000 -s 4 eval -o data/results_4.json
    cargo run -r -- -n 1000000 -s 256 eval -o data/results_256.json

fast:
    # cargo run -r -- -n 500000 -s 4 eval -o data/results_4.json
    cargo run -r -- -n 500000 -s 256 eval -o data/results_256.json

run:
    cargo run -r -- -n 100000 -s 256 eval

stat:
    cargo build -r
    perf stat cargo run -r -- -n 100000 -s 256 eval

flame:
    cargo flamegraph --open --skip-after minimizers::MinimizerType::density -- -n 100000 -s 256 eval
