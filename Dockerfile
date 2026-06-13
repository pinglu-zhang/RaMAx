FROM ubuntu:22.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_JOBS=4

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        cmake \
        g++-12 \
        gcc-12 \
        libcurl4-openssl-dev \
        libhdf5-dev \
        libtbb-dev \
        make \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

ENV CC=/usr/bin/gcc-12
ENV CXX=/usr/bin/g++-12

WORKDIR /src

COPY CMakeLists.txt ./
COPY include ./include
COPY src ./src
COPY third_party ./third_party

RUN set -eux; \
    hdf5_cflags="$(pkg-config --cflags hdf5)"; \
    hdf5_libs="$(pkg-config --libs hdf5)"; \
    export CFLAGS="${hdf5_cflags}"; \
    export CXXFLAGS="${hdf5_cflags}"; \
    cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/opt/ramax \
        -DCMAKE_INSTALL_LIBDIR=/opt/ramax/lib \
        -DCMAKE_INSTALL_PKGCONFIGDIR=/opt/ramax/lib/pkgconfig \
        -DCMAKE_INSTALL_RPATH=/opt/ramax/lib \
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
        -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
        -DBUILD_EXAMPLES=OFF \
        -DBUILD_SHARED_LIBS=OFF \
        -DRAMAX_NATIVE_ARCH=OFF \
        -DRAMAX_HAL_JOBS="${BUILD_JOBS}" \
        -DRAMAX_HAL_LIBS="${hdf5_libs} -lhdf5_cpp"; \
    cmake --build build --parallel "${BUILD_JOBS}"; \
    cmake --install build; \
    /opt/ramax/bin/ramax --help >/dev/null; \
    if /opt/ramax/bin/ramax --help | grep -q -- "--mask-repeats"; then exit 1; fi; \
    test ! -e /src/bin; \
    test ! -e /opt/ramax/bin/windowmasker

FROM ubuntu:22.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        libcurl4 \
        libgomp1 \
        libhdf5-103-1 \
        libhdf5-cpp-103-1 \
        libstdc++6 \
        libtbb12 \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/ramax /opt/ramax

RUN set -eux; \
    /opt/ramax/bin/ramax --help >/dev/null; \
    if /opt/ramax/bin/ramax --help | grep -q -- "--mask-repeats"; then exit 1; fi; \
    ldd /opt/ramax/bin/ramax; \
    if ldd /opt/ramax/bin/ramax | grep -q "not found"; then exit 1; fi; \
    test ! -e /opt/ramax/bin/windowmasker; \
    ! command -v windowmasker

WORKDIR /data

ENTRYPOINT ["/opt/ramax/bin/ramax"]
CMD ["--help"]
