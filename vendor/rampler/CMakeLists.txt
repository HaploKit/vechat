cmake_minimum_required(VERSION 3.9)

project(rampler VERSION 2.0.1
                LANGUAGES CXX
                DESCRIPTION "Rampler is a tool for subsampling or splitting FASTA/Q files.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (NOT TARGET bioparser)
  add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
endif ()
if (NOT TARGET biosoup)
  add_subdirectory(vendor/bioparser/vendor/biosoup EXCLUDE_FROM_ALL)
endif ()

add_executable(${PROJECT_NAME}
  src/main.cpp
  src/sampler.cpp)
target_link_libraries(${PROJECT_NAME}
  bioparser
  biosoup)
target_compile_definitions(${PROJECT_NAME} PRIVATE
  RAMPLER_VERSION="v${PROJECT_VERSION}")

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_BINDIR})
