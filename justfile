eval:
    cargo run -r -- -n 300000 -s 4 eval > results_4.json
    cargo run -r -- -n 300000 -s 256 eval > results_256.json

run:
    cargo run -r -- -n 10000 -s 4 eval
