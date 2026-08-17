cmake_minimum_required(VERSION 3.20)

project(RaMAx VERSION 1.0.5 LANGUAGES C CXX)
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

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/ramax_version.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/ramax_version.h"
    @ONLY)

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
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Local builds can opt into native CPU tuning for better performance. Conda
# packages disable this in recipe/build.sh so the binary remains portable.
option(RAMAX_NATIVE_ARCH "Build ramax with -march=native" ON)

set(RAMAX_MINIPOA_EXECUTABLE "" CACHE FILEPATH
    "Existing minipoa executable; empty searches PATH")

if(RAMAX_MINIPOA_EXECUTABLE)
    if(NOT EXISTS "${RAMAX_MINIPOA_EXECUTABLE}")
        message(FATAL_ERROR
            "RAMAX_MINIPOA_EXECUTABLE does not exist: "
            "${RAMAX_MINIPOA_EXECUTABLE}")
    endif()
    get_filename_component(RAMAX_RESOLVED_MINIPOA
        "${RAMAX_MINIPOA_EXECUTABLE}" REALPATH)
    message(STATUS "Using configured minipoa: ${RAMAX_RESOLVED_MINIPOA}")
else()
    find_program(RAMAX_SYSTEM_MINIPOA NAMES minipoa)
    if(RAMAX_SYSTEM_MINIPOA)
        set(RAMAX_RESOLVED_MINIPOA "${RAMAX_SYSTEM_MINIPOA}")
        message(STATUS "Using minipoa from PATH: ${RAMAX_RESOLVED_MINIPOA}")
    else()
        set(RAMAX_RESOLVED_MINIPOA "")
        message(STATUS
            "minipoa was not found during configuration; install it in PATH "
            "or set RAMAX_MINIPOA_EXECUTABLE before using MSA-based options")
    endif()
endif()

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build bundled libraries statically" FORCE)
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory("${THIRD_PARTY_DIR}/libdivsufsort" EXCLUDE_FROM_ALL)

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
    "${ALIGN_DIR}/cross_anchor_repair.cpp"
    "${ALIGN_DIR}/external_msa_runner.cpp"
    "${ALIGN_DIR}/reference_profile_merger.cpp"
    "${ALIGN_DIR}/cigar.cpp"
    "${ALIGN_DIR}/ksw2_extz2_sse.c"
    "${ALIGN_DIR}/pair_rare_aligner.cpp"
    "${ALIGN_DIR}/multiple_rare_aligner.cpp"
    "${ALIGN_DIR}/minipoa_locator.cpp"
    "${GRAPH_DIR}/ramesh_node.cpp"
    "${GRAPH_DIR}/ramesh_graph.cpp"
    "${GRAPH_DIR}/direct_block_merger.cpp"
    "${GRAPH_DIR}/realignment_internal.cpp"
    "${GRAPH_DIR}/structural_break_repair.cpp"
    "${GRAPH_DIR}/short_block_repair.cpp"
    "${GRAPH_DIR}/ramesh_export.cpp"
    "${GRAPH_DIR}/paf_utils.cpp"
    "${GRAPH_DIR}/paf_export.cpp"
    "${GRAPH_DIR}/hal_converter.cpp"
)

set(MAIN_X "${SOURCE_DIR}/RaMAx.cpp")

# ================= 公共头文件路径 =================
set(RAMA_INCLUDE_DIRS
    "${INCLUDE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/generated"
)

set(RAMA_SYSTEM_INCLUDE_DIRS
    "${THIRD_PARTY_DIR}/spdlog/include"
    "${BIGBWT_INCLUDE_DIR}"
    "${SUFFIX_ARRAY_INCLUDE_DIR}"
    "${CMAKE_BINARY_DIR}/third_party/libdivsufsort/include"
    "${THIRD_PARTY_DIR}/parlaylib/include"
    "${HAL_DIR}/api/inc"
    "${SONLIB_DIR}/C/inc"
    ${HDF5_CXX_INCLUDE_DIRS}
)

# ================= 编译配置函数 =================
function(configure_rama_target target)
    target_include_directories(${target} PRIVATE ${RAMA_INCLUDE_DIRS})
    target_include_directories(${target} SYSTEM PRIVATE ${RAMA_SYSTEM_INCLUDE_DIRS})

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

    add_dependencies(${target} divsufsort divsufsort64)

    target_compile_definitions(${target} PRIVATE M64)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_definitions(${target} PRIVATE _DEBUG_)
    endif()

    target_compile_options(${target} PRIVATE
        -Wall
        -Wextra
        -Wno-unused-parameter
        # kseq and bigBWT are embedded upstream headers included by legacy code.
        -Wno-unused-function
        -Wno-sign-compare
        # GCC 13 emits a false positive for Match vector copy elision.
        -Wno-array-bounds
        -Wno-stringop-overflow
        -Wno-deprecated-declarations
        -Wno-builtin-macro-redefined)
    if(RAMAX_NATIVE_ARCH)
        target_compile_options(${target} PRIVATE -march=native)
    endif()

endfunction()

# Vendored algorithm implementations retain their upstream warning policy.
set_source_files_properties(
    "${INDEX_DIR}/bigBWT_build.cpp"
    "${INDEX_DIR}/bigBWT/newscan.cpp"
    "${INDEX_DIR}/bigBWT/utils.c"
    "${INDEX_DIR}/bigBWT/xerrors.c"
    "${SUFFIX_ARRAY_DIR}/gsa/gsacak.c"
    "${SUFFIX_ARRAY_DIR}/CaPS-SA/Suffix_Array.cpp"
    "${ANCHOR_DIR}/xxhash.c"
    "${ALIGN_DIR}/ksw2_extz2_sse.c"
    PROPERTIES COMPILE_OPTIONS "-w")

# ================= shared RaMAx core =================
add_library(ramax_core STATIC ${RAMA_COMMON_SOURCES})
configure_rama_target(ramax_core)
add_dependencies(ramax_core hal_library)
target_compile_definitions(ramax_core PRIVATE
    RAMAX_MINIPOA_CONFIGURED_PATH="${RAMAX_RESOLVED_MINIPOA}")

add_executable(ramax ${MAIN_X})
configure_rama_target(ramax)
target_link_libraries(ramax PRIVATE ramax_core)

add_dependencies(ramax  hal_library)

include(CTest)
if(BUILD_TESTING)
    add_executable(ramax_output_spec_tests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/output_spec_test.cpp")
    target_include_directories(ramax_output_spec_tests PRIVATE "${INCLUDE_DIR}")
    target_compile_features(ramax_output_spec_tests PRIVATE cxx_std_23)
    target_compile_options(ramax_output_spec_tests PRIVATE -Wall -Wextra)
    add_test(NAME ramax_output_spec COMMAND ramax_output_spec_tests)

    add_test(NAME ramax_multi_output_cli
        COMMAND "${CMAKE_COMMAND}"
            -DRAMAX=$<TARGET_FILE:ramax>
            -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/multi-output-cli-test
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/multi_output_cli_test.cmake")

    add_executable(ramax_paf_tests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/paf_utils_test.cpp"
        "${GRAPH_DIR}/paf_utils.cpp")
    target_include_directories(ramax_paf_tests PRIVATE "${INCLUDE_DIR}")
    target_compile_features(ramax_paf_tests PRIVATE cxx_std_23)
    target_compile_options(ramax_paf_tests PRIVATE -Wall -Wextra)
    add_test(NAME ramax_paf_utils COMMAND ramax_paf_tests)

    add_executable(ramax_paf_export_tests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/paf_export_test.cpp")
    configure_rama_target(ramax_paf_export_tests)
    target_link_libraries(ramax_paf_export_tests PRIVATE ramax_core)
    target_compile_features(ramax_paf_export_tests PRIVATE cxx_std_23)
    add_test(NAME ramax_paf_export COMMAND ramax_paf_export_tests)
endif()

install(TARGETS ramax DESTINATION "${CMAKE_INSTALL_BINDIR}")

# Optional companion tools are self-contained under tools/.  A tool that was
# already built in its own source-local build directory is installed directly,
# avoiding a second compilation as part of the main RaMAx build.
option(RAMAX_BUILD_TOOLS "Build and install RaMAx companion tools" ON)
set(RAMAX_PAF_FASTA_PREBUILT "" CACHE FILEPATH
    "Prebuilt ramax-paf-fasta executable; skips building the tool")

if(RAMAX_BUILD_TOOLS)
    set(_ramax_paf_fasta_prebuilt "${RAMAX_PAF_FASTA_PREBUILT}")
    if(NOT _ramax_paf_fasta_prebuilt)
        set(_ramax_paf_fasta_auto
            "${CMAKE_CURRENT_SOURCE_DIR}/tools/ramax-paf-fasta/build/ramax-paf-fasta${CMAKE_EXECUTABLE_SUFFIX}")
        if(EXISTS "${_ramax_paf_fasta_auto}")
            set(_ramax_paf_fasta_prebuilt "${_ramax_paf_fasta_auto}")
        endif()
    elseif(NOT EXISTS "${_ramax_paf_fasta_prebuilt}")
        message(FATAL_ERROR
            "RAMAX_PAF_FASTA_PREBUILT does not exist: ${_ramax_paf_fasta_prebuilt}")
    endif()

    if(_ramax_paf_fasta_prebuilt)
        message(STATUS
            "Using prebuilt ramax-paf-fasta: ${_ramax_paf_fasta_prebuilt}; source build skipped")
        install(PROGRAMS "${_ramax_paf_fasta_prebuilt}"
            DESTINATION "${CMAKE_INSTALL_BINDIR}"
            RENAME "ramax-paf-fasta${CMAKE_EXECUTABLE_SUFFIX}")
    else()
        message(STATUS "ramax-paf-fasta is not prebuilt; adding standalone tool build")
        add_subdirectory(tools/ramax-paf-fasta)
    endif()
endif()
