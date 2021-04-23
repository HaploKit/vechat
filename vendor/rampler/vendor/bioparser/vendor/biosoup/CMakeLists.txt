cmake_minimum_required(VERSION 3.9)

project(biosoup VERSION 0.10.0
                LANGUAGES CXX
                DESCRIPTION "Biosoup is a c++ collection of header only data structures used for storage and logging in bioinformatics tools.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

option(biosoup_build_tests "Build biosoup unit tests" OFF)
if (biosoup_build_tests)
  find_package(GTest REQUIRED)
  add_executable(${PROJECT_NAME}_test
    test/nucleic_acid_test.cpp
    test/overlap_test.cpp
    test/progress_bar_test.cpp
    test/sequence_test.cpp
    test/timer_test.cpp)
  target_link_libraries(${PROJECT_NAME}_test
    ${PROJECT_NAME}
    GTest::Main)
endif ()
