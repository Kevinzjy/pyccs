#!/bin/bash

yum -y update
yum -y groupinstall "Development Tools"
yum -y install wget git openssl-devel openssh-client 
    
cd /tmp
wget https://github.com/Kitware/CMake/releases/download/v3.20.1/cmake-3.20.1.tar.gz
tar zxvf cmake-3.20.1.tar.gz
cd cmake-3.20.1
./configure
gmake
make install

cd /tmp 
wget https://zlib.net/fossils/zlib-1.2.11.tar.gz
tar zxvf zlib-1.2.11.tar.gz
cd zlib-1.2.11
./configure
make 
make install

export PATH=/opt/cross/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH
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
