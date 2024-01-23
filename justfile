eval:
    cargo run -r -- -n 100000 -s 4 eval > results_4.json
    cargo run -r -- -n 100000 -s 16 eval > results_16.json
    cargo run -r -- -n 100000 -s 256 eval > results_256.json
