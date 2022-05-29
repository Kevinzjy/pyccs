#!/bin/bash
source /root/.bashrc

export PATH=/opt/cross/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/include:$LD_LIBRARY_PATH
export ZLIB_LIBRARY=/usr/local/lib
export ZLIB_INCLUDE_DIR=/usr/local/include

cp ./misc/config /root/.cargo/config

cargo clean

# Build wheelss
for PYBIN in /opt/python/cp3[7]*/bin; do
    "${PYBIN}/pip" install maturin -i https://pypi.tuna.tsinghua.edu.cn/simple
    "${PYBIN}/maturin" build -i "${PYBIN}/python" --release
done

# Repair wheels
for wheel in target/wheels/*.whl; do
    auditwheel repair "${wheel}"
done
