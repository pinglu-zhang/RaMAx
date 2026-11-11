#!/usr/bin/env bash
set -euo pipefail

export CFLAGS="${CFLAGS:-} -O3 -I${PREFIX}/include"
export CXXFLAGS="${CXXFLAGS:-} -O3 -I${PREFIX}/include"
export CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS:-}"
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PREFIX}/share/pkgconfig:${PKG_CONFIG_PATH:-}"
test -e "${PREFIX}/lib/libhdf5.so"
test -e "${PREFIX}/lib/libhdf5_cpp.so"
hal_libs="${LDFLAGS:-} -Wl,-rpath,${PREFIX}/lib -Wl,-rpath-link,${PREFIX}/lib -L${PREFIX}/lib -lhdf5_cpp -lhdf5"

cmake -S . -B build-conda \
  ${CMAKE_ARGS:-} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DCMAKE_INSTALL_LIBDIR="${PREFIX}/lib" \
  -DCMAKE_INSTALL_PKGCONFIGDIR="${PREFIX}/lib/pkgconfig" \
  -DCMAKE_INSTALL_RPATH="${PREFIX}/lib" \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF -DRAMAX_BUILD_TESTS=OFF \
  -DRAMAX_INSTALL_VALIDATION_TOOLS=ON \
  -DRAMAX_BUILD_TOOLS=ON -DRAMAX_EMBED_TOOL_PATHS=OFF -DRAMAX_NATIVE_ARCH=OFF \
  -DRAMAX_HAL_LIBS="${hal_libs}" -DRAMAX_HAL_JOBS="${CPU_COUNT:-2}"
cmake --build build-conda --parallel "${CPU_COUNT:-2}"
cmake --install build-conda

test "$("${PREFIX}/bin/ramax" --version)" = "RaMAx version ${PKG_VERSION}"
test -x "${PREFIX}/bin/ramax-paf-fasta"
test -x "${PREFIX}/bin/ramax-spawn-helper"
test -x "${PREFIX}/libexec/ramax/ramax_release_hal_check"
install -d "${PREFIX}/share/licenses/ramax" "${PREFIX}/share/ramax"
install -m 0644 LICENSE "${PREFIX}/share/licenses/ramax/RaMAx-LICENSE"
install -m 0644 third_party/hal/LICENSE.txt "${PREFIX}/share/licenses/ramax/HAL-LICENSE.txt"
{
  printf 'tool\tversion_or_source\n'
  printf 'ramax\t%s\n' "${PKG_VERSION}"
  printf 'HAL\tnative bundled library\n'
  printf 'mash\t2.3=hb105d93_9\n'
  printf 'wfmash\t0.14.0=h11f254b_0\n'
  printf 'samtools\t1.23.1=ha83d96e_0\n'
  printf 'htslib\t1.23.1=h633afcb_0\n'
  printf 'minipoa\t1.4.2=hd5d28ae_0\n'
} > "${PREFIX}/share/ramax/packaged-tools.tsv"
