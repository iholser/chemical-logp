# docker build -t ianh/chem-notebook .
# Takes ~2.5hrs to build

FROM debian:bullseye-slim

RUN apt-get update && apt-get install -y \
  apt-utils \
  build-essential \
  cmake \
  git \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-system-dev \
  libboost-thread-dev \
  libboost-serialization-dev \
  libboost-python-dev \
  libboost-regex-dev \
  libeigen3-dev \
  libfreetype6 \
  libfreetype6-dev \
  libpng-dev \
  libsqlite3-dev \
  pkg-config \
  python3-numpy \
  python3-dev \
  python3-pip \
  python3-matplotlib \
  sqlite3 \
  wget \
  zip \
  && apt-get upgrade -y

RUN pip3 install --upgrade pip setuptools

ENV RDKIT_BRANCH=Release_2021_09
RUN git clone -b $RDKIT_BRANCH --single-branch https://github.com/rdkit/rdkit.git

ENV RDBASE=/rdkit
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RDBASE/lib:/usr/lib/x86_64-linux-gnu
ENV PYTHONPATH=$PYTHONPATH:$RDBASE

RUN mkdir $RDBASE/build
WORKDIR $RDBASE/build

RUN cmake \
  -D RDK_BUILD_INCHI_SUPPORT=ON \
  ..
RUN make && make install

ADD ./requirements.txt /tmp/requirements.txt
RUN pip3 install -qr /tmp/requirements.txt

WORKDIR $RDBASE
