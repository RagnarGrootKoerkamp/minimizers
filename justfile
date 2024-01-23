eval-small:
    cargo run -r -- -n 10000 -s 4 eval > results_4.json
    # cargo run -r -- -n 10000 -s 16 eval > results_16.json
    # cargo run -r -- -n 10000 -s 256 eval > results_256.json
eval:
    cargo run -r -- -n 100000 -s 4 eval > results_4.json
    # cargo run -r -- -n 100000 -s 16 eval > results_16.json
    # cargo run -r -- -n 100000 -s 256 eval > results_256.json

run:
    cargo run -r -- -n 100000 -s 4 eval
