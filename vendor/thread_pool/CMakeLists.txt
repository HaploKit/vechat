cmake_minimum_required(VERSION 3.9)

project(thread_pool VERSION 3.0.2
                    LANGUAGES CXX
                    DESCRIPTION "ThreadPool is a c++ header only library modifying and extending https://github.com/progschj/ThreadPool.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

find_package(Threads REQUIRED)
add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE
  Threads::Threads)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

option(thread_pool_build_tests "Build thread_pool unit tests" OFF)
if (thread_pool_build_tests)
  find_package(GTest REQUIRED)
  add_executable(${PROJECT_NAME}_test
    test/semaphore_test.cpp
    test/thread_pool_test.cpp)
  target_link_libraries(${PROJECT_NAME}_test
    ${PROJECT_NAME}
    GTest::Main)
endif ()
