cmake_minimum_required(VERSION 3.20)

project(ramax-paf-fasta VERSION 1.0.5 LANGUAGES CXX)

include(GNUInstallDirs)
find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)

add_library(ramax_paf_fasta_lib STATIC
    src/ramax_paf_fasta.cpp)
target_include_directories(ramax_paf_fasta_lib
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(ramax_paf_fasta_lib PUBLIC cxx_std_23)
target_compile_options(ramax_paf_fasta_lib PRIVATE -Wall -Wextra)
target_link_libraries(ramax_paf_fasta_lib PUBLIC CURL::libcurl ZLIB::ZLIB)

add_executable(ramax-paf-fasta src/main.cpp)
target_compile_features(ramax-paf-fasta PRIVATE cxx_std_23)
target_compile_options(ramax-paf-fasta PRIVATE -Wall -Wextra)
target_compile_definitions(ramax-paf-fasta PRIVATE
    RAMAX_PAF_FASTA_VERSION="${PROJECT_VERSION}")
target_link_libraries(ramax-paf-fasta PRIVATE ramax_paf_fasta_lib)

install(TARGETS ramax-paf-fasta DESTINATION "${CMAKE_INSTALL_BINDIR}")
