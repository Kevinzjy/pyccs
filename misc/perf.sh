#!/bin/bash
# RUSTFLAGS=-g cargo build --release
cargo build --release
# perf record ./target/release/ccs -t 8 -i ./tests/test.fasta.gz -o ./tests/ccs.fa -r ./tests/raw.fa
time ./target/release/ccs -t 8 -i ./tests/sim_circ.fasta.gz -o ./tests/ccs.fa -r ./tests/raw.fa
