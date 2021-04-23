# Bioparser

[![Latest GitHub release](https://img.shields.io/github/release/rvaser/bioparser.svg)](https://github.com/rvaser/bioparser/releases/latest)
[![Build status for gcc/clang](https://travis-ci.com/rvaser/bioparser.svg?branch=master)](https://travis-ci.com/rvaser/bioparser)

Bioparser is a c++ header only parsing library for several formats in bioinformatics (FASTA/Q, MHAP/PAF/SAM), with support for zlib compressed files.

## Usage

If you would like to add bioparser to your project via CMake, add the following:
```cmake
if (NOT TARGET bioparser)
  add_subdirectory(<path_to_submodules>/bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(<your_exe> bioparser)
```

If you are not using CMake, include the appropriate header file directly to your project and link with zlib.

#### Dependencies
- gcc 4.8+ or clang 3.5+
- (optional) cmake 3.9+
- zlib

## Examples

### FASTA parser

```cpp
#include "bioparser/fasta_parser.hpp"

struct Sequence {  // or any other name
 public:
  Sequence(  // required arguments
      const char*, std::uint32_t,
      const char*, std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(path);

// parse whole file
auto s = p->Parse(-1);
```

### FASTQ parser

```cpp
#include "bioparser/fastq_parser.hpp"

struct Sequence {  // or any other name
 public:
  Sequence(  // required arguments
      const char*, std::uint32_t,
      const char*, std::uint32_t,
      const char*, std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(path);

// parse in chunks
std::vector<std::unique_ptr<Sequence>> s;
std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
for (auto t = p->parse(chunk_size); !t.empty(); t = p->parse(chunk_size)) {
  s.insert(
      s.end(),
      std::make_move_iterator(t.begin()),
      std::make_move_iterator(t.end()));
}
```

### MHAP parser

```cpp
#include "bioparser/mhap_parser.hpp"

struct Overlap {  // or any other name
 public:
  Overlap(  // required arguments
      std::uint64_t,
      std::uint64_t,
      double error,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Overlap>::Create<bioparser::MhapParser>(path);

// parse whole file
auto o = p->Parse(-1);
```

### PAF parser

```cpp
#include "bioparser/paf_parser.hpp"

struct Overlap {  // or any other name
 public:
  Overlap(  // required arguments
      const char*, std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      char,
      const char*, std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Overlap>::Create<bioparser::PafParser>(path);

// parse whole file
auto o = p->Parse(-1);
```

### SAM parser

```cpp
#include "bioparser/sam_parser.hpp"

struct Overlap {  // or any other name
 public:
  Overlap(  // required arguments
      const char*, std::uint32_t,
      std::uint32_t,
      const char*, std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      const char*, std::uint32_t,
      const char*, std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      const char*, std::uint32_t,
      const char*, std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Overlap>::Create<bioparser::SamParser>(path);

// parse whole file
auto o = p->Parse(-1);
```

**Note**: If your class has a private constructor, add one of the following lines to your class definition:

```cpp
friend bioparser::FastaParser<Sequence>;
friend bioparser::FastqParser<Sequence>;
friend bioparser::MhapParser<Overlap>;
friend bioparser::PafParser<Overlap>;
friend bioparser::SamParser<Overlap>;
```

## Unit tests

To build and run bioparser unit tests run the following commands:

```bash
git clone --recursive https://github.com/rvaser/bioparser.git bioparser
cd bioparser && mkdir build && cd build
cmake -Dbioparser_build_tests=ON -DCMAKE_BUILD_TYPE=Release .. && make
./bin/bioparser_test
```

#### Dependencies
- gtest

## Acknowledgement

This work has been supported in part by the Croatian Science Foundation under the project Single genome and metagenome assembly (IP-2018-01-5886).
