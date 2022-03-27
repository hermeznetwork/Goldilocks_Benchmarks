FROM ubuntu:21.10

WORKDIR /usr/src/app

RUN apt update && apt install -y \
    build-essential \
    libgmp-dev \
    libbenchmark-dev \
    nasm \
    libomp-dev

COPY . /usr/src/app

RUN make -j

ENTRYPOINT [ "/usr/src/app/build/benchmark" ]
