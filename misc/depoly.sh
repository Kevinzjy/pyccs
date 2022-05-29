#!/bin/bash

# git clone --recursive https://github.com/Kevinzjy/pyccs.git
# cd pyccs/misc

# podman build -t manylinux2010_x86_64_rust:v1 ./

# setproxy

# podman run -v \
#   /home/zhangjy/data/git/pyccs:/volume:z \
#   -w /volume \
#   --userns keep-id \
#   --network host \
#   --rm \
#   -t localhost/centos7_rust:v1 cargo test --verbose && cargo build --verbose --release

# strip ./target/release/ccs

# mv ./target/release/ccs ./target/release/ccs_v1.0.0_el7.x86-64

# RUSTFLAGS=-g cargo build --release
# cargo build --release
# perf record ./target/release/ccs -t 8 -i ./test/test.fq.gz -o ./tests/ccs.fa -r ./tests/raw.fa

setproxy

podman run -v \
  /home/zhangjy/data/git/pyccs:/volume:z \
  -w /volume \
  --network host \
  --rm \
  -i \
  -t localhost/manylinux2010_x86_64_rust:v1 
  sh ./misc/build-wheels.sh 


  -i \

  --userns keep-id \