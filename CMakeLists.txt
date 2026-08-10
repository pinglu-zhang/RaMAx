cmake_minimum_required(VERSION 3.10)

project("RaMAx")
if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif()

# MSVC 热重载支持
if(POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT
    "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,"
    "$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

# 设置 C++23 标准
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ================= 路径设置 =================
set(SOURCE_DIR         "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INCLUDE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(THIRD_PARTY_DIR    "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
set(PREPROCESS_DIR     "${SOURCE_DIR}/preprocess")
set(INDEX_DIR          "${SOURCE_DIR}/index")
set(ALIGN_DIR          "${SOURCE_DIR}/align")
set(ANCHOR_DIR         "${SOURCE_DIR}/anchor")
set(GRAPH_DIR          "${SOURCE_DIR}/graph")
set(BIGBWT_DIR         "${INDEX_DIR}/bigBWT")
set(BIGBWT_INCLUDE_DIR "${INCLUDE_DIR}/bigBWT")
set(SUFFIX_ARRAY_DIR   "${INDEX_DIR}/suffix_array")
set(SUFFIX_ARRAY_INCLUDE_DIR "${INCLUDE_DIR}/suffix_array")
set(SONLIB_DIR ${THIRD_PARTY_DIR}/sonLib)
set(HAL_DIR    ${THIRD_PARTY_DIR}/hal)

find_program(MAKE_EXE NAMES make gmake REQUIRED)
set(RAMAX_HAL_JOBS "" CACHE STRING "Parallel jobs for bundled sonLib/HAL make builds")
set(RAMAX_HAL_LIBS "" CACHE STRING "Additional linker flags for bundled HAL tools")
if(RAMAX_HAL_JOBS)
    set(RAMAX_HAL_MAKE_JOBS "-j${RAMAX_HAL_JOBS}")
else()
    set(RAMAX_HAL_MAKE_JOBS "-j")
endif()

# ================= 外部库依赖 =================
find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)

set(RAMAX_HAL_HDF5_INCLUDE_FLAGS "")
foreach(dir IN LISTS HDF5_CXX_INCLUDE_DIRS)
    set(RAMAX_HAL_HDF5_INCLUDE_FLAGS "${RAMAX_HAL_HDF5_INCLUDE_FLAGS} -I${dir}")
endforeach()

# ================= HAL library =================
add_custom_command(
    OUTPUT ${SONLIB_DIR}/lib/sonLib.h ${SONLIB_DIR}/lib/sonLib.a
    COMMAND ${MAKE_EXE} -C ${SONLIB_DIR} ${RAMAX_HAL_MAKE_JOBS} libs
            "CC=${CMAKE_C_COMPILER}"
            "CXX=${CMAKE_CXX_COMPILER}"
            "AR=${CMAKE_AR}"
            "RANLIB=${CMAKE_RANLIB}"
    WORKING_DIRECTORY ${SONLIB_DIR}
    COMMENT "Building sonLib…"
)

add_custom_target(sonLib_build
    DEPENDS ${SONLIB_DIR}/lib/sonLib.a
)

add_custom_command(
    OUTPUT ${HAL_DIR}/lib/libHal.a
    COMMAND ${MAKE_EXE} -C ${HAL_DIR} ${RAMAX_HAL_MAKE_JOBS} libs
            "CC=${CMAKE_C_COMPILER}"
            "CXX=${CMAKE_CXX_COMPILER}"
            "AR=${CMAKE_AR}"
            "RANLIB=${CMAKE_RANLIB}"
            "HAL_HDF5_CFLAGS=${RAMAX_HAL_HDF5_INCLUDE_FLAGS}"
            "LIBS=${RAMAX_HAL_LIBS}"
    WORKING_DIRECTORY ${HAL_DIR}
    COMMENT "Building HAL library..."
)

add_custom_target(hal_library
    DEPENDS ${HAL_DIR}/lib/libHal.a
)

add_dependencies(hal_library sonLib_build)
# 如果没有显式设置构建类型，默认使用 Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Local builds can opt into native CPU tuning for better performance. Conda
# packages disable this in recipe/build.sh so the binary remains portable.
option(RAMAX_NATIVE_ARCH "Build ramax with -march=native" ON)

# 去掉符号表（减小二进制体积）
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-s")


# ================= bundled dependencies =================
# add_subdirectory("${THIRD_PARTY_DIR}/CLI11")
# add_subdirectory("${THIRD_PARTY_DIR}/sdsl-lite")
#set(WFA_BUILD_SHARED OFF CACHE BOOL "" FORCE)
#set(WFA_BUILD_STATIC ON CACHE BOOL "" FORCE)
# add_subdirectory("${THIRD_PARTY_DIR}/WFA2-lib")
# 已经 add_subdirectory(...) 产生了目标 "divsufsort"
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build bundled libraries statically" FORCE)
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory("${THIRD_PARTY_DIR}/libdivsufsort")

# ================= 公共源文件 =================
set(RAMA_COMMON_SOURCES
    "${PREPROCESS_DIR}/raw_data_process.cpp"
    "${PREPROCESS_DIR}/newick_tree.cpp"
    "${PREPROCESS_DIR}/sequence_utils.cpp"
    "${PREPROCESS_DIR}/SeqPro.cpp"
    "${PREPROCESS_DIR}/softmask_index.cpp"
    "${INDEX_DIR}/fm_index.cpp"
    "${INDEX_DIR}/bigBWT_build.cpp"
    "${INDEX_DIR}/bigBWT/newscan.cpp"
    "${INDEX_DIR}/bigBWT/utils.c"
    "${INDEX_DIR}/bigBWT/xerrors.c"
    "${SUFFIX_ARRAY_DIR}/gsa/gsacak.c"
    "${SUFFIX_ARRAY_DIR}/CaPS-SA/Suffix_Array.cpp"
    "${ANCHOR_DIR}/anchor.cpp"
    "${ANCHOR_DIR}/anchor_chunk.cpp"
    "${ANCHOR_DIR}/cluster.cpp"
    "${ANCHOR_DIR}/xxhash.c"
    "${ALIGN_DIR}/align.cpp"
    "${ALIGN_DIR}/cigar.cpp"
    "${ALIGN_DIR}/ksw2_extz2_sse.c"
    "${ALIGN_DIR}/pair_rare_aligner.cpp"
    "${ALIGN_DIR}/multiple_rare_aligner.cpp"
    "${GRAPH_DIR}/ramesh_node.cpp"
    "${GRAPH_DIR}/ramesh_graph.cpp"
    "${GRAPH_DIR}/direct_block_merger.cpp"
    "${GRAPH_DIR}/window_detector.cpp"
    "${GRAPH_DIR}/ramesh_export.cpp"
    "${GRAPH_DIR}/hal_converter.cpp"
)

# ================= 主程序入口文件 =================
#set(MAIN_G "${SOURCE_DIR}/RaMA-G.cpp")
set(MAIN_X "${SOURCE_DIR}/RaMAx.cpp")

# ================= 公共头文件路径 =================
set(RAMA_INCLUDE_DIRS
    "${INCLUDE_DIR}"
    "${THIRD_PARTY_DIR}/spdlog/include"
    # "${THIRD_PARTY_DIR}/CLI11/include"
#   "${THIRD_PARTY_DIR}/cereal/include"
    "${BIGBWT_INCLUDE_DIR}"
    "${SUFFIX_ARRAY_INCLUDE_DIR}"
    # "${THIRD_PARTY_DIR}/sdsl-lite/include"
    # "${THIRD_PARTY_DIR}/sdsl-lite/external/libdivsufsort/include"
    "${CMAKE_BINARY_DIR}/third_party/libdivsufsort/include"
    # "${THIRD_PARTY_DIR}/WFA2-lib"
    "${THIRD_PARTY_DIR}/parlaylib/include"
    "${HAL_DIR}/api/inc"
    "${SONLIB_DIR}/C/inc"
    ${HDF5_CXX_INCLUDE_DIRS}
)

# ================= 编译配置函数 =================
function(configure_rama_target target)
    target_include_directories(${target} PRIVATE ${RAMA_INCLUDE_DIRS})

    target_link_libraries(${target} PRIVATE
        # CLI11::CLI11
        ZLIB::ZLIB
        CURL::libcurl
        # sdsl
        divsufsort
        divsufsort64
        Threads::Threads
        #wfa2cpp_static
        TBB::tbb
        OpenMP::OpenMP_CXX
        ${HAL_DIR}/lib/libHal.a
        ${SONLIB_DIR}/lib/sonLib.a
        ${HDF5_CXX_LIBRARIES}
    )

    # 确保在构建目标前先构建libdivsufsort（生成头文件）
    add_dependencies(${target} divsufsort divsufsort64)

    target_compile_definitions(${target} PRIVATE M64 SPDLOG_HEADER_ONLY)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_definitions(${target} PRIVATE _DEBUG_)
    endif()

    target_compile_options(${target} PRIVATE -Wall -g)
    if(RAMAX_NATIVE_ARCH)
        target_compile_options(${target} PRIVATE -march=native)
    endif()

endfunction()

# add_compile_definitions(KSW_CPU_DISPATCH)
add_compile_options(-w)

# ================= 创建两个可执行文件 =================
# add_executable(RaMA-G ${MAIN_G} ${RAMA_COMMON_SOURCES})
# clonfigure_rama_target(RaMA-G)

add_executable(ramax ${MAIN_X} ${RAMA_COMMON_SOURCES})
configure_rama_target(ramax)

# 如需确保先有 HAL 工具，再编你的程序：
#add_dependencies(RaMA-G hal_tools)
add_dependencies(ramax  hal_library)

# ================= 安装规则（可选） =================
# install(TARGETS RaMAG RaMAx DESTINATION bin)

# ================= 安装规则 =================
# 安装主程序
install(TARGETS ramax DESTINATION bin)

# ================= optional regression tests =================
option(RAMAX_BUILD_TESTS "Build RaMAx soft-mask regression tests" OFF)
if(RAMAX_BUILD_TESTS)
    enable_testing()

    add_executable(ramax_softmask_index_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/softmask_index_test.cpp"
        "${PREPROCESS_DIR}/softmask_index.cpp")
    target_include_directories(ramax_softmask_index_test PRIVATE
        "${INCLUDE_DIR}")
    target_link_libraries(ramax_softmask_index_test PRIVATE ZLIB::ZLIB)
    target_compile_options(ramax_softmask_index_test PRIVATE -Wall)
    add_test(NAME ramax_softmask_index COMMAND ramax_softmask_index_test)

    add_executable(ramax_hal_softmask_roundtrip_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/hal_softmask_roundtrip_test.cpp"
        "${PREPROCESS_DIR}/softmask_index.cpp")
    configure_rama_target(ramax_hal_softmask_roundtrip_test)
    add_dependencies(ramax_hal_softmask_roundtrip_test hal_library)
    add_test(NAME ramax_hal_softmask_roundtrip
        COMMAND ramax_hal_softmask_roundtrip_test)

    add_executable(ramax_hal_softmask_export_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/hal_softmask_export_test.cpp"
        ${RAMA_COMMON_SOURCES})
    configure_rama_target(ramax_hal_softmask_export_test)
    add_dependencies(ramax_hal_softmask_export_test hal_library)
    add_test(NAME ramax_hal_softmask_export
        COMMAND ramax_hal_softmask_export_test)

    add_executable(ramax_window_detector_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/window_detector_test.cpp"
        ${RAMA_COMMON_SOURCES})
    configure_rama_target(ramax_window_detector_test)
    add_dependencies(ramax_window_detector_test hal_library)
    add_test(NAME ramax_window_detector COMMAND ramax_window_detector_test)

    add_executable(ramax_zero_gap_missing_species_merge_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/zero_gap_missing_species_merge_test.cpp"
        ${RAMA_COMMON_SOURCES})
    configure_rama_target(ramax_zero_gap_missing_species_merge_test)
    add_dependencies(ramax_zero_gap_missing_species_merge_test hal_library)
    add_test(NAME ramax_zero_gap_missing_species_merge
        COMMAND ramax_zero_gap_missing_species_merge_test)
endif()
