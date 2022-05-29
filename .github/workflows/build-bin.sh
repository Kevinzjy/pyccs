#!/bin/bash   
export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/include:$LD_LIBRARY_PATH
export ZLIB_LIBRARY=/usr/local/lib
export ZLIB_INCLUDE_DIR=/usr/local/include

PYBIN=/opt/python/cp37-cp37m/bin
${PYBIN}/maturin build -i ${PYBIN}/python --release -b bin
