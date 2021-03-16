#!/bin/bash
# RUSTFLAGS=-g cargo build --release
cargo build --release
time ./target/release/ccs -t 8 -i ./tests/simulated_reads.fasta.gz -o ./tests/ccs.fa -r ./tests/raw.fa
