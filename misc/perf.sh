#!/bin/bash
# RUSTFLAGS=-g cargo build --release
cargo build --release
perf record ./target/release/ccs -t 8 -i ./test/test.fq.gz -o ./tests/ccs.fa -r ./tests/raw.fa
