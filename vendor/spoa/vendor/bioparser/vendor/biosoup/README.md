# Biosoup

[![Build status for gcc/clang](https://travis-ci.com/rvaser/biosoup.svg?branch=master)](https://travis-ci.com/rvaser/biosoup)

Biosoup is a c++ collection of header only data structures used for storage and logging in bioinformatics tools.

## Usage

If you would like to add biosoup to your project via CMake, add the following:
```cmake
if (NOT TARGET biosoup)
  add_subdirectory(<path_to_submodules>/biosoup EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(<your_exe> biosoup)
```

If you are not using CMake, include the appropriate header file directly to your project.

#### Dependencies

- gcc 4.8+ or clang 3.5+
- (optional) cmake 3.9+

## Unit tests

To build and run biosoup unit tests run the following commands:

```bash
git clone https://github.com/rvaser/biosoup.git biosoup
cd biosoup && mkdir build && cd build
cmake -Dbiosoup_build_tests=ON -DCMAKE_BUILD_TYPE=Release .. && make
./bin/biosoup_test
```

#### Dependencies
- gtest

## Acknowledgement

This work has been supported in part by the Croatian Science Foundation under the project Single genome and metagenome assembly (IP-2018-01-5886).
