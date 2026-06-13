#!/usr/bin/env bash
set -euxo pipefail

export CFLAGS="${CFLAGS:-} -O3 -I${PREFIX}/include"
export CXXFLAGS="${CXXFLAGS:-} -O3 -I${PREFIX}/include"
export CPPFLAGS="-I${PREFIX}/include ${CPPFLAGS:-}"
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PREFIX}/share/pkgconfig:${PKG_CONFIG_PATH:-}"

# The source checkout may contain HAL/sonLib artifacts built outside conda.
# Reusing them can cause HDF5 C++ ABI/link mismatches, so force conda to
# rebuild these static libraries and tools inside the isolated build env.
rm -rf \
  build-conda \
  third_party/hal/bin \
  third_party/hal/lib \
  third_party/hal/objs \
  third_party/sonLib/bin \
  third_party/sonLib/lib \
  third_party/sonLib/C/impl/*.o \
  third_party/sonLib/externalTools/quicktree_1.1/obj

# HAL normally compiles through h5cc/h5c++, which can bypass conda-build's
# compiler/sysroot selection. The top-level CMake build passes CC/CXX through
# to HAL and uses these flags to retain HDF5 and conda-prefix linkage.
HDF5_LINK_FLAGS="$(pkg-config --libs hdf5)"
RAMAX_HAL_LIBS="${LDFLAGS:-} -Wl,-rpath,${PREFIX}/lib -Wl,-rpath-link,${PREFIX}/lib -L${PREFIX}/lib -lhdf5_cpp ${HDF5_LINK_FLAGS}"

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
  -DBUILD_SHARED_LIBS=OFF \
  -DRAMAX_NATIVE_ARCH=OFF \
  -DRAMAX_HAL_LIBS="${RAMAX_HAL_LIBS}" \
  -DRAMAX_HAL_JOBS="${RAMAX_HAL_JOBS:-4}"

cmake --build build-conda --parallel "${CPU_COUNT:-2}"
cmake --install build-conda
