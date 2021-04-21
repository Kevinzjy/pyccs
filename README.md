# circtools

Accelerating functions in CIRI toolkit.

For any questions, please contact zhangjinyang@biols.ac.cn

## Release Notes

- Version 1.0: First released version

## License

The code is released under the `MIT` License. See the LICENSE file for more detail

## Installation 

### Binary distribution

Compiled binary can be downloaded from [Github Release](https://github.com/Kevinzjy/circtools/releases)

Download the binary according to your platform. For example:

```bash
wget -O ccs https://github.com/Kevinzjy/circtools/releases/download/v1.0.0/ccs_v1.0.0_el7.x86_64
chmod +x ccs
./ccs --help
```

### Install from source 

#### Dependency

- rust 1.50.0+
- cmake 3.12+
- zlib 1.2.8+

#### Compile circtools

1. Install Rust 

```bash
curl https://sh.rustup.rs -sSf | sh
```

2. Install circtools

```bash
git clone --recursive https://github.com/Kevinzjy/circtools.git
cargo build --release
```

3. Copy compiled binary to somewhere you want, and add it to `$PATH`

```bash
cp ./target/release/ccs your/path
export PATH=your/path:$PATH
```

## Q&A

**(1) Building CMake Error: Could NOT find ZLIB**

Install zlib 1.2.8+, and add the path of your zlib library in `build.rs`.

```rust
let dst = Config::new("vendor/spoa")
        .define("CMAKE_BUILD_TYPE","Release")
        .define("ZLIB_LIBRARY", "/your/path/to/zlib/lib")
        .define("ZLIB_INCLUDE_DIR", "/your/path/to/zlib/include")
        .build();
```

Then, run `cargo build --release` to build `ccs` with your custom zlib version.

**(2) /lib64/libc.so.6: version `GLIBC_2.18' not found not found**

For example, when running on machines with `GLIBC<2.18`:

```
./ccs_v1.0.0_x86_64-unknown-linux-gnu_generic: /lib64/libc.so.6: version `GLIBC_2.18' not found (required by ./ccs_v1.0.0_x86_64-unknown-linux-gnu_generic)
./ccs_v1.0.0_x86_64-unknown-linux-gnu_generic: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by ./ccs_v1.0.0_x86_64-unknown-linux-gnu_generic)
./ccs_v1.0.0_x86_64-unknown-linux-gnu_generic: /lib64/libstdc++.so.6: version `CXXABI_1.3.8' not found (required by ./ccs_v1.0.0_x86_64-unknown-linux-gnu_generic)
```

The released `ccs` binaries are compiled using Ubuntu 16.04, which use `GLIBC_2.18` as the default glibc version. If you're running ccs on a linux machine with older GLIBC library, please compile the `ccs` binary from source. Please refer to the instructions above.
