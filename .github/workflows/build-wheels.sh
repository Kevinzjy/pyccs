#!/bin/bash   
export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/include:$LD_LIBRARY_PATH
export ZLIB_LIBRARY=/usr/local/lib
export ZLIB_INCLUDE_DIR=/usr/local/include

for PYBIN in /opt/python/cp3[7891]*/bin; do
    "${PYBIN}/pip" install maturin
    "${PYBIN}/maturin" build -i "${PYBIN}/python" --release
done

for wheel in target/wheels/*.whl; do
    auditwheel repair "${wheel}"
done
