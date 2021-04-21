#!/bin/bash
# Install on CentOS 7 with glibc 2.17

# Install dependencies
strings /lib64/libc.so.6 | grep -P "^GLIBC"
yum groupinstall "Development Tools"
yum install git openssl-devel

# Update cmake
wget https://github.com/Kitware/CMake/releases/download/v3.20.1/cmake-3.20.1.tar.gz
tar zxvf cmake-3.20.1.tar.gz
cd cmake-3.20.1
configure && gmake && make install 

wget http://zlib.net/zlib-1.2.11.tar.gz
tar zxvf zlib-1.1.11.tar.gz
cd zlib-1.2.11
configure && make && make install

# Install Rust
useradd rust
su - rust
curl https://sh.rustup.rs -sSf | sh

# Compile circtools
git clone --recursive https://github.com/Kevinzjy/circtools.git
cd circtools
cargo build --release
mv ./target/release/ccs ./target/release/ccs_v1.0.0_el7.x86-64

# RUSTFLAGS=-g cargo build --release
# cargo build --release
# perf record ./target/release/ccs -t 8 -i ./test/test.fq.gz -o ./tests/ccs.fa -r ./tests/raw.fa
