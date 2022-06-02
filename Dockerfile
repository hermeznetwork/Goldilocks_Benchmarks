FROM ubuntu:21.10

WORKDIR /usr/src/app

RUN apt update && apt install -y build-essential libomp-dev cmake git libgmp3-dev

RUN git clone https://github.com/google/benchmark.git && \
    cd benchmark  && \
    cmake -E make_directory "build"   && \
    cmake -E chdir "build" cmake -DBENCHMARK_DOWNLOAD_DEPENDENCIES=on -DCMAKE_BUILD_TYPE=Release ../ && \
    cmake --build "build" --config Release && \ 
    cmake --build "build" --config Release --target install

COPY . /usr/src/app

RUN make -j

ENTRYPOINT [ "/usr/src/app/build/benchmark" ]
