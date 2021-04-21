FROM centos:7
RUN yum -y update && \
    yum -y groupinstall "Development Tools" && \
    yum -y install wget git openssl-devel 
RUN mkdir /rust
WORKDIR /rust
RUN wget https://github.com/Kitware/CMake/releases/download/v3.20.1/cmake-3.20.1.tar.gz \
    && tar zxvf cmake-3.20.1.tar.gz \
    && cd cmake-3.20.1 \
    && ./configure \
    && gmake \
    && make install
RUN wget http://zlib.net/zlib-1.2.11.tar.gz \
    && tar zxvf zlib-1.2.11.tar.gz \
    && cd zlib-1.2.11 \
    && ./configure \
    && make \
    && make install

RUN useradd rust --user-group --create-home --shell /bin/bash --groups wheel
RUN echo "rust ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
USER rust
RUN mkdir -p /home/rust/libs /home/rust/src
ENV PATH=/home/rust/.cargo/bin:/opt/cross/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
RUN sh -c "curl https://sh.rustup.rs -sSf | sh -s -- -y"
RUN git config --global credential.https://github.com.helper ghtoken
WORKDIR /home/rust/src
