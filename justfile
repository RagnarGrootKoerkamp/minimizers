eval:
    cargo run -r -- -n 5000000 -s 256 eval -o data/density_256.json
eval-4:
    cargo run -r -- -n 5000000 -s 4 eval -o data/density_4.json

stats:
    cargo run -r -- -n 5000000 -s 256 eval --stats -o data/stats_256.json

flame:
    cargo flamegraph --open --skip-after minimizers::MinimizerType::stats -- -n 100000 -s 256 eval
flame-stat:
    cargo flamegraph --open --skip-after minimizers::MinimizerType::stats -- -n 100000 -s 256 eval --stats

py:
    source .env/bin/activate && maturin develop -r
