# Rampler

[![Latest GitHub release](https://img.shields.io/github/release/rvaser/rampler.svg)](https://github.com/rvaser/rampler/releases/latest)
[![Build status for c++/clang++](https://travis-ci.com/rvaser/rampler.svg?branch=master)](https://travis-ci.com/rvaser/rampler)

Rampler is a standalone module for sampling genomic sequences. It supports two modes, random subsampling of sequencing data to a desired depth (given the reference length) and file splitting to desired size in bytes.

## Usage

To build rampler run the following commands:
```bash
git clone --recursive https://github.com/rvaser/rampler.git rampler
cd rampler && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
./bin/rampler
```

which will display the following usage:

```bash
usage: rampler [options ...] <mode>

  <mode>
    subsample <sequences> <reference length> <coverage> [<coverage> ...]

      <sequences>
        input file in FASTA/FASTQ format (can be compressed with gzip)
      <reference length>
        integer denoting length of the reference genome (or assembly)
      <coverage>
        integer denoting desired coverage of the subsampled sequences

    split <sequences> <chunk size>

      <sequences>
        input file in FASTA/FASTQ format (can be compressed with gzip)
        containing sequences which will be split into smaller chunks
      <chunk size>
        integer denoting the desired chunk size in bytes

  options:
    -o, --out-directory <string>
      default: current directory
      path in which sampled files will be created
    --version
      prints the version number
    -h, --help
      prints out the help
```

#### Dependencies
1. gcc 4.8+ or clang 4.0+
2. cmake 3.9+
3. zlib

## Acknowledgment
This work has been supported in part by Croatian Science Foundation under the project UIP-11-2013-7353.
