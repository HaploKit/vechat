cmake_minimum_required(VERSION 3.9)

project(bioparser VERSION 3.0.12
                  LANGUAGES CXX
                  DESCRIPTION "Bioparser is a c++ header only parsing library for several formats in bioinformatics (FASTA/Q, MHAP/PAF/SAM), with support for zlib compressed files.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

find_package(ZLIB REQUIRED)
add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE
  ZLIB::ZLIB)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

option(bioparser_build_tests "Build bioparser unit tests" OFF)
if (bioparser_build_tests)
  find_package(GTest REQUIRED)
  if (NOT TARGET biosoup)
    add_subdirectory(vendor/biosoup EXCLUDE_FROM_ALL)
  endif ()
  add_executable(${PROJECT_NAME}_test
    test/parser_test.cpp
    test/fasta_parser_test.cpp
    test/fastq_parser_test.cpp
    test/mhap_parser_test.cpp
    test/paf_parser_test.cpp
    test/sam_parser_test.cpp)
  target_link_libraries(${PROJECT_NAME}_test
    ${PROJECT_NAME}
    biosoup
    GTest::Main)
  target_compile_definitions(${PROJECT_NAME}_test PRIVATE
    BIOPARSER_DATA_PATH="${PROJECT_SOURCE_DIR}/test/data/")
endif ()
