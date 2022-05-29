#!/bin/bash
set -e -x

source /root/.bashrc

# Clone repository
# cd /tmp
# git clone git@github.com:Kevinzjy/pyccs.git
# cd pyccs

# Build wheelss
for PYBIN in /opt/python/cp3[6]*/bin; do
    "${PYBIN}/pip" install maturin
    "${PYBIN}/maturin" build -i "${PYBIN}/python" --release
done

# Repair wheels
for wheel in target/wheels/*.whl; do
    auditwheel repair "${wheel}"
done

ls /wheelhouse/*.whl
# mv /wheelhouse/*.whl /volume/target/wheels

#  ENV PATH=/home/rust/.cargo/bin:/opt/cross/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin