# circtools

Accelerating functions in CIRI toolkit

## Author

Author & Maintainer: Jinyang Zhang(zhangjinyang@biols.ac.cn)

## Release Notes

version 1.0: First released version

## License

The code is released under the MIT License. See the LICENSE file for more detail

## Installation 

### Dependency

- rust 1.50.0
- cmake 3.2+

### Compile circtools

1. Install Rust 

```bash
curl https://sh.rustup.rs -sSf | sh
```

2. Install circtools

```bash
git clone --recursive https://github.com/Kevinzjy/circtools.git
cargo build --release
```

The process should be finished within ~1min.

3. Copy compiled binary to somewhere you want, and add it to $PATH

```bash
cp ./target/release/ccs your/path
export PATH=your/path:$PATH
```

