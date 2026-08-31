cmake_minimum_required(VERSION 3.20)

project(RaMAx VERSION 1.0.9 LANGUAGES C CXX)
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
option(RAMAX_EMBED_TOOL_PATHS
    "Embed resolved external-tool paths in the RaMAx binary" ON)
set(RAMAX_TOOL_BIN_DIR "${CMAKE_SOURCE_DIR}/bin" CACHE PATH
    "Directory containing RaMAx external tools for source builds")
if(NOT RAMAX_TOOL_BIN_DIR OR NOT IS_DIRECTORY "${RAMAX_TOOL_BIN_DIR}")
    message(FATAL_ERROR
        "RAMAX_TOOL_BIN_DIR is not an existing directory: ${RAMAX_TOOL_BIN_DIR}")
endif()
get_filename_component(RAMAX_RESOLVED_TOOL_BIN_DIR
    "${RAMAX_TOOL_BIN_DIR}" REALPATH)
message(STATUS
    "RaMAx external-tool directory: ${RAMAX_RESOLVED_TOOL_BIN_DIR}")

set(RAMAX_MINIPOA_EXECUTABLE "" CACHE FILEPATH
    "Explicit minipoa executable; empty uses the shared tool lookup order")
set(RAMAX_MASH_EXECUTABLE "" CACHE FILEPATH
    "Explicit Mash 2.3 executable; empty uses the shared tool lookup order")
set(RAMAX_SAMTOOLS_EXECUTABLE "" CACHE FILEPATH
    "Explicit Samtools 1.23.1 executable; empty uses shared tool lookup")
set(RAMAX_WFMASH_EXECUTABLE "" CACHE FILEPATH
    "Explicit PGGB-compatible wfmash v0.14 executable; empty uses shared lookup")
set(RAMAX_HAL_APPEND_CACTUS_SUBTREE_EXECUTABLE "" CACHE FILEPATH
    "Explicit halAppendCactusSubtree executable; empty uses shared lookup")

function(resolve_ramax_tool configured_variable tool_name display_name output_variable)
    if(${configured_variable})
        if(NOT EXISTS "${${configured_variable}}")
            message(FATAL_ERROR
                "${configured_variable} does not exist: ${${configured_variable}}")
        endif()
        get_filename_component(_resolved "${${configured_variable}}" REALPATH)
        message(STATUS "Using configured ${display_name}: ${_resolved}")
    else()
        set(_resolved "")
        set(_bundled_candidate
            "${RAMAX_RESOLVED_TOOL_BIN_DIR}/${tool_name}")
        if(EXISTS "${_bundled_candidate}")
            message(STATUS
                "Using ${display_name} from RAMAX_TOOL_BIN_DIR: "
                "${_bundled_candidate}")
        else()
            find_program(_system_tool NAMES "${tool_name}" NO_CACHE)
            if(_system_tool)
                message(WARNING
                    "${display_name} is absent from RAMAX_TOOL_BIN_DIR; it is "
                    "currently available through PATH at ${_system_tool}, but "
                    "that PATH result will not be embedded")
            else()
                message(WARNING
                    "${display_name} is absent from RAMAX_TOOL_BIN_DIR and "
                    "PATH; it must be provided before the corresponding "
                    "runtime stage")
            endif()
        endif()
    endif()
    set(${output_variable} "${_resolved}" PARENT_SCOPE)
endfunction()

resolve_ramax_tool(RAMAX_MASH_EXECUTABLE mash "Mash 2.3"
                   RAMAX_RESOLVED_MASH)
resolve_ramax_tool(RAMAX_SAMTOOLS_EXECUTABLE samtools "Samtools 1.23.1"
                   RAMAX_RESOLVED_SAMTOOLS)
resolve_ramax_tool(RAMAX_WFMASH_EXECUTABLE wfmash "wfmash v0.14.0-0-g517e1bc"
                   RAMAX_RESOLVED_WFMASH)
resolve_ramax_tool(RAMAX_HAL_APPEND_CACTUS_SUBTREE_EXECUTABLE
                   halAppendCactusSubtree "halAppendCactusSubtree"
                   RAMAX_RESOLVED_HAL_APPEND_CACTUS_SUBTREE)

resolve_ramax_tool(RAMAX_MINIPOA_EXECUTABLE minipoa "minipoa"
                   RAMAX_RESOLVED_MINIPOA)

set(RAMAX_COMPILED_MINIPOA "${RAMAX_RESOLVED_MINIPOA}")
set(RAMAX_COMPILED_MASH "${RAMAX_RESOLVED_MASH}")
set(RAMAX_COMPILED_SAMTOOLS "${RAMAX_RESOLVED_SAMTOOLS}")
set(RAMAX_COMPILED_WFMASH "${RAMAX_RESOLVED_WFMASH}")
set(RAMAX_COMPILED_HAL_APPEND_CACTUS_SUBTREE
    "${RAMAX_RESOLVED_HAL_APPEND_CACTUS_SUBTREE}")
set(RAMAX_COMPILED_TOOL_BIN_DIR "${RAMAX_RESOLVED_TOOL_BIN_DIR}")
if(NOT RAMAX_EMBED_TOOL_PATHS)
    set(RAMAX_COMPILED_MINIPOA "")
    set(RAMAX_COMPILED_MASH "")
    set(RAMAX_COMPILED_SAMTOOLS "")
    set(RAMAX_COMPILED_WFMASH "")
    set(RAMAX_COMPILED_HAL_APPEND_CACTUS_SUBTREE "")
    set(RAMAX_COMPILED_TOOL_BIN_DIR "")
    message(STATUS
        "External-tool paths will not be embedded; runtime lookup uses the "
        "ramax executable directory and PATH")
endif()

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build bundled libraries statically" FORCE)
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory("${THIRD_PARTY_DIR}/libdivsufsort" EXCLUDE_FROM_ALL)

# ================= 公共源文件 =================
set(RAMA_COMMON_SOURCES
    "${PREPROCESS_DIR}/cache_manifest.cpp"
    "${PREPROCESS_DIR}/raw_data_process.cpp"
    "${PREPROCESS_DIR}/newick_tree.cpp"
    "${PREPROCESS_DIR}/sequence_utils.cpp"
    "${PREPROCESS_DIR}/SeqPro.cpp"
    "${PREPROCESS_DIR}/softmask_index.cpp"
    "${INDEX_DIR}/fm_index.cpp"
    "${INDEX_DIR}/suffix_array_index.cpp"
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
    "${ALIGN_DIR}/external_tool.cpp"
    "${ALIGN_DIR}/dependency_preflight.cpp"
    "${ALIGN_DIR}/reference_profile_merger.cpp"
    "${ALIGN_DIR}/cigar.cpp"
    "${ALIGN_DIR}/ksw2_extz2_sse.c"
    "${ALIGN_DIR}/pair_rare_aligner.cpp"
    "${ALIGN_DIR}/multiple_rare_aligner.cpp"
    "${ALIGN_DIR}/minipoa_locator.cpp"
    "${ALIGN_DIR}/mash_distance_estimator.cpp"
    "${ALIGN_DIR}/wfmash_router.cpp"
    "${GRAPH_DIR}/ramesh_node.cpp"
    "${GRAPH_DIR}/ramesh_graph.cpp"
    "${GRAPH_DIR}/direct_block_merger.cpp"
    "${GRAPH_DIR}/realignment_internal.cpp"
    "${GRAPH_DIR}/structural_break_repair.cpp"
    "${GRAPH_DIR}/short_block_repair.cpp"
    "${GRAPH_DIR}/ramesh_export.cpp"
    "${GRAPH_DIR}/aligned_block_view.cpp"
    "${GRAPH_DIR}/gfa_export.cpp"
    "${SOURCE_DIR}/hal/export.cpp"
    "${GRAPH_DIR}/paf_utils.cpp"
    "${GRAPH_DIR}/paf_export.cpp"
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
    RAMAX_TOOL_BIN_CONFIGURED_PATH="${RAMAX_COMPILED_TOOL_BIN_DIR}"
    RAMAX_MINIPOA_CONFIGURED_PATH="${RAMAX_COMPILED_MINIPOA}"
    RAMAX_MASH_CONFIGURED_PATH="${RAMAX_COMPILED_MASH}"
    RAMAX_SAMTOOLS_CONFIGURED_PATH="${RAMAX_COMPILED_SAMTOOLS}"
    RAMAX_WFMASH_CONFIGURED_PATH="${RAMAX_COMPILED_WFMASH}"
    RAMAX_HAL_APPEND_CACTUS_SUBTREE_CONFIGURED_PATH="${RAMAX_COMPILED_HAL_APPEND_CACTUS_SUBTREE}")

add_executable(ramax ${MAIN_X})
configure_rama_target(ramax)
target_link_libraries(ramax PRIVATE ramax_core)

add_dependencies(ramax hal_library)

option(BUILD_TESTING "Build CTest-based RaMAx tests" OFF)
include(CTest)

if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    add_executable(ramax_mash_distance_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/mash_distance_estimator_test.cpp")
    configure_rama_target(ramax_mash_distance_test)
    target_link_libraries(ramax_mash_distance_test PRIVATE ramax_core)
    add_test(NAME ramax_mash_distance_test
        COMMAND ramax_mash_distance_test)
    set_tests_properties(ramax_mash_distance_test PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

    add_executable(ramax_wfmash_router_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/wfmash_router_test.cpp")
    configure_rama_target(ramax_wfmash_router_test)
    target_link_libraries(ramax_wfmash_router_test PRIVATE ramax_core)
    add_test(NAME ramax_wfmash_router_test
        COMMAND ramax_wfmash_router_test)
    set_tests_properties(ramax_wfmash_router_test PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

    add_executable(ramax_gfa_export_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/gfa_export_test.cpp")
    configure_rama_target(ramax_gfa_export_test)
    target_link_libraries(ramax_gfa_export_test PRIVATE ramax_core)
    add_test(NAME ramax_gfa_export_test
        COMMAND ramax_gfa_export_test)
    set_tests_properties(ramax_gfa_export_test PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/core_graph_performance_equivalence_test.cpp")
        add_executable(ramax_core_graph_performance_equivalence_test
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/core_graph_performance_equivalence_test.cpp")
        configure_rama_target(ramax_core_graph_performance_equivalence_test)
        target_compile_definitions(ramax_core_graph_performance_equivalence_test
            PRIVATE RAMAX_PERFORMANCE_TEST_HOOKS)
        target_link_libraries(ramax_core_graph_performance_equivalence_test
            PRIVATE ramax_core)
        add_test(NAME ramax_core_graph_performance_equivalence_test
            COMMAND ramax_core_graph_performance_equivalence_test)
        set_tests_properties(ramax_core_graph_performance_equivalence_test PROPERTIES
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    endif()
endif()

install(TARGETS ramax DESTINATION "${CMAKE_INSTALL_BINDIR}")

option(RAMAX_BUILD_TESTS "Build RaMAx regression tests" OFF)
if(RAMAX_BUILD_TESTS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    add_executable(ramax_hal_export_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/hal_export_test.cpp")
    configure_rama_target(ramax_hal_export_test)
    target_link_libraries(ramax_hal_export_test PRIVATE ramax_core)
    add_dependencies(ramax_hal_export_test hal_library)
    add_test(NAME ramax_hal_export COMMAND ramax_hal_export_test)

    add_executable(ramax_hal_softmask_export_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/hal_softmask_export_test.cpp")
    configure_rama_target(ramax_hal_softmask_export_test)
    target_link_libraries(ramax_hal_softmask_export_test PRIVATE ramax_core)
    add_dependencies(ramax_hal_softmask_export_test hal_library)
    add_test(NAME ramax_hal_softmask_export
        COMMAND ramax_hal_softmask_export_test)

    add_executable(ramax_zero_gap_missing_species_merge_test
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/zero_gap_missing_species_merge_test.cpp"
        ${RAMA_COMMON_SOURCES})
    configure_rama_target(ramax_zero_gap_missing_species_merge_test)
    target_compile_definitions(ramax_zero_gap_missing_species_merge_test
        PRIVATE RAMAX_KSW_TEST_HOOKS)
    add_dependencies(ramax_zero_gap_missing_species_merge_test hal_library)
    add_test(NAME ramax_zero_gap_missing_species_merge
        COMMAND ramax_zero_gap_missing_species_merge_test)
endif()

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
