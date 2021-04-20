#!/bin/bash
# RUSTFLAGS=-g cargo build --release
cargo build --release
# perf record ./target/release/ccs

# Create test files
gunzip -c test.fa.gz > test.fa
ln -s test.fa test.fasta
ln -s test.fa.gz test.fasta.gz
gunzip -c test.fq.gz > test.fq
ln -s test.fq test.fastq
ln -s test.fq.gz test.fastq.gz

for suffix in "fa" "fasta" "fa.gz" "fasta.gz" "fq" "fastq" "fq.gz" "fastq.gz"
do
    ../target/release/ccs -t 8 -i test.${suffix} -o ccs.fa -r raw.fa
done

rm ccs.fa raw.fa
rm test.fa test.fq
rm test.fasta test.fasta.gz
rm test.fastq test.fastq.gz
