ARG RAMAX_VERSION=1.0.9

FROM mambaorg/micromamba:2.3.2 AS alignment-tools

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root
RUN micromamba create -y --strict-channel-priority \
        -p /opt/ramax \
        -c conda-forge -c bioconda -c malab \
        mash=2.3=hb105d93_9 \
        minipoa=1.4.2=hd5d28ae_0 \
        wfmash=0.14.0=h11f254b_0 \
        samtools=1.23.1=ha83d96e_0 \
        htslib=1.23.1=h633afcb_0 \
    && micromamba clean --all --yes \
    && test "$(/opt/ramax/bin/mash --version)" = "2.3" \
    && test "$(/opt/ramax/bin/minipoa -v)" = "minipoa version 1.4.2" \
    && test "$(/opt/ramax/bin/wfmash --version 2>&1 | sed -n '1p')" = "v0.14.0-0-g517e1bc" \
    && test "$(/opt/ramax/bin/samtools --version | sed -n '1p')" = "samtools 1.23.1" \
    && test "$(/opt/ramax/bin/samtools --version | sed -n '2p')" = "Using htslib 1.23.1"

FROM ubuntu:22.04 AS cactus-helper

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        binutils \
        file \
        libc-bin \
        patchelf \
    && rm -rf /var/lib/apt/lists/*

RUN install -d /opt/ramax/bin /opt/ramax/share/ramax
COPY bin/halAppendCactusSubtree \
    /opt/ramax/bin/halAppendCactusSubtree

RUN set -eux; \
    helper=/opt/ramax/bin/halAppendCactusSubtree; \
    chmod 0755 "${helper}"; \
    description="$(file "${helper}")"; \
    printf '%s\n' "${description}"; \
    printf '%s\n' "${description}" | grep -q 'ELF 64-bit.*x86-64'; \
    readelf -d "${helper}"; \
    objdump -T "${helper}" | grep -o 'GLIBC_[0-9.]*' | sort -Vu; \
    if printf '%s\n' "${description}" | grep -q 'statically linked'; then \
        helper_help="$("${helper}" --help 2>&1 || true)"; \
        printf '%s\n' "${helper_help}" | grep -q 'USAGE:'; \
    else \
        before_ldd="$(ldd "${helper}")"; \
        printf '%s\n' "${before_ldd}"; \
        if printf '%s\n' "${before_ldd}" | grep -qE 'not found|/mnt/sda/'; then exit 1; fi; \
        if patchelf --print-rpath "${helper}" | grep -q '/mnt/sda/'; then patchelf --remove-rpath "${helper}"; fi; \
        after_ldd="$(ldd "${helper}")"; \
        printf '%s\n' "${after_ldd}"; \
        if printf '%s\n' "${after_ldd}" | grep -qE 'not found|/mnt/sda/'; then exit 1; fi; \
    fi

FROM ubuntu:22.04 AS builder

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_JOBS=4
ARG RAMAX_VERSION

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        cmake \
        g++-12 \
        gcc-12 \
        libboost-graph-dev \
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

COPY --from=alignment-tools /opt/ramax /opt/ramax
COPY --from=cactus-helper /opt/ramax /opt/ramax

COPY CMakeLists.txt ./
COPY LICENSE THIRD_PARTY_NOTICES.md ./
COPY include ./include
COPY src ./src
COPY tools ./tools
COPY third_party ./third_party

RUN set -eux; \
    test ! -e /src/tests; \
    test ! -e /src/tools/ramax-paf-fasta/build/ramax-paf-fasta; \
    test "$(sed -n 's/^project(RaMAx VERSION \([^ ]*\).*/\1/p' CMakeLists.txt)" = "${RAMAX_VERSION}"; \
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
        -DBUILD_TESTING=OFF \
        -DRAMAX_BUILD_TESTS=OFF \
        -DRAMAX_BUILD_TOOLS=ON \
        -DRAMAX_TOOL_BIN_DIR=/opt/ramax/bin \
        -DRAMAX_EMBED_TOOL_PATHS=ON \
        -DRAMAX_NATIVE_ARCH=OFF \
        -DRAMAX_MASH_EXECUTABLE=/opt/ramax/bin/mash \
        -DRAMAX_WFMASH_EXECUTABLE=/opt/ramax/bin/wfmash \
        -DRAMAX_SAMTOOLS_EXECUTABLE=/opt/ramax/bin/samtools \
        -DRAMAX_MINIPOA_EXECUTABLE=/opt/ramax/bin/minipoa \
        -DRAMAX_HAL_APPEND_CACTUS_SUBTREE_EXECUTABLE=/opt/ramax/bin/halAppendCactusSubtree \
        -DRAMAX_HAL_JOBS="${BUILD_JOBS}" \
        -DRAMAX_HAL_LIBS="${hdf5_libs} -lhdf5_cpp"; \
    cmake --build build --parallel "${BUILD_JOBS}"; \
    cmake --install build; \
    install -d /opt/ramax/share/licenses/ramax /opt/ramax/share/ramax; \
    install -m 0644 LICENSE /opt/ramax/share/licenses/ramax/RaMAx-LICENSE; \
    install -m 0644 third_party/hal/LICENSE.txt /opt/ramax/share/licenses/ramax/HAL-LICENSE; \
    /opt/ramax/bin/ramax --help >/dev/null; \
    test "$(/opt/ramax/bin/ramax --version)" = "RaMAx version ${RAMAX_VERSION}"; \
    test -x /opt/ramax/bin/ramax-paf-fasta; \
    if /opt/ramax/bin/ramax --help | grep -q -- "--mask-repeats"; then exit 1; fi; \
    test ! -e /src/bin; \
    test ! -e /opt/ramax/bin/windowmasker

FROM ubuntu:22.04 AS runtime

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
ARG RAMAX_VERSION
ARG VCS_REF=unknown

LABEL org.opencontainers.image.title="RaMAx" \
      org.opencontainers.image.description="Whole-genome alignment and pangenome graph construction" \
      org.opencontainers.image.version="${RAMAX_VERSION}" \
      org.opencontainers.image.source="https://github.com/pinglu-zhang/RaMAx" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.vendor="malab"

ENV PATH=/opt/ramax/bin:${PATH}

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
    test "$(/opt/ramax/bin/ramax --version)" = "RaMAx version ${RAMAX_VERSION}"; \
    /opt/ramax/bin/ramax --help >/dev/null; \
    if /opt/ramax/bin/ramax --help | grep -q -- "--mask-repeats"; then exit 1; fi; \
    ldd /opt/ramax/bin/ramax | tee /tmp/ramax.ldd; \
    if grep -q "not found" /tmp/ramax.ldd; then exit 1; fi; \
    for executable in ramax mash wfmash samtools minipoa halAppendCactusSubtree; do \
        test "$(command -v "${executable}")" = "/opt/ramax/bin/${executable}"; \
    done; \
    helper_ldd="$(ldd /opt/ramax/bin/halAppendCactusSubtree 2>&1 || true)"; \
    printf '%s\n' "${helper_ldd}"; \
    if printf '%s\n' "${helper_ldd}" | grep -q 'not a dynamic executable'; then \
        helper_help="$(halAppendCactusSubtree --help 2>&1 || true)"; \
        printf '%s\n' "${helper_help}" | grep -q 'USAGE:'; \
    elif printf '%s\n' "${helper_ldd}" | grep -qE 'not found|/mnt/sda/'; then exit 1; fi; \
    test "$(mash --version)" = "2.3"; \
    test "$(wfmash --version 2>&1 | sed -n '1p')" = "v0.14.0-0-g517e1bc"; \
    test "$(samtools --version | sed -n '1p')" = "samtools 1.23.1"; \
    test "$(samtools --version | sed -n '2p')" = "Using htslib 1.23.1"; \
    test "$(minipoa -v)" = "minipoa version 1.4.2"; \
    test ! -e /opt/ramax/bin/windowmasker; \
    test ! -e /opt/ramax/tests; \
    ! command -v windowmasker

WORKDIR /data

ENTRYPOINT ["/opt/ramax/bin/ramax"]
CMD ["--help"]
