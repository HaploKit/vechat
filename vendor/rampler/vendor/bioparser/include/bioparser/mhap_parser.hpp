// Copyright (c) 2020 Robert Vaser

#ifndef BIOPARSER_MHAP_PARSER_HPP_
#define BIOPARSER_MHAP_PARSER_HPP_

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <vector>
#include <stdexcept>

#include "bioparser/parser.hpp"

namespace bioparser {

template<class T>
class MhapParser: public Parser<T> {
 public:
  MhapParser(const MhapParser&) = delete;
  MhapParser& operator=(const MhapParser&) = delete;

  MhapParser(MhapParser&&) = delete;
  MhapParser& operator=(MhapParser&&) = delete;

  ~MhapParser() {}

  std::vector<std::unique_ptr<T>> Parse(
      std::uint64_t bytes, bool = true) override {
    std::vector<std::unique_ptr<T>> dst;
    std::uint64_t parsed_bytes = 0;

    std::uint64_t lhs_id = 0;
    std::uint64_t rhs_id = 0;
    double error = 0;
    std::uint32_t num_minmers = 0;
    std::uint32_t lhs_strand = 0;
    std::uint32_t lhs_begin = 0;
    std::uint32_t lhs_end = 0;
    std::uint32_t lhs_len = 0;
    std::uint32_t rhs_strand = 0;
    std::uint32_t rhs_begin = 0;
    std::uint32_t rhs_end = 0;
    std::uint32_t rhs_len = 0;

    auto create_T = [&] () -> void {
      auto storage_ptr = this->RightStrip(
          this->storage().data(),
          this->storage_ptr());
      this->Terminate(storage_ptr);

      std::uint32_t num_values = 0;
      std::uint32_t begin_ptr = 0;
      while (true) {
        auto end_ptr = begin_ptr;
        while (end_ptr < storage_ptr && this->storage()[end_ptr] != ' ') {
          ++end_ptr;
        }
        this->Terminate(end_ptr);

        switch (num_values) {
          case 0: lhs_id = std::atoll(this->storage().data() + begin_ptr); break;  // NOLINT
          case 1: rhs_id = std::atoll(this->storage().data() + begin_ptr); break;  // NOLINT
          case 2: error = std::atof(this->storage().data() + begin_ptr); break;
          case 3: num_minmers = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 4: lhs_strand = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 5: lhs_begin = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 6: lhs_end = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 7: lhs_len = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 8: rhs_strand = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 9: rhs_begin = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 10: rhs_end = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 11: rhs_len = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          default: break;
        }

        ++num_values;
        if (end_ptr == storage_ptr || num_values == 12) {
          break;
        }
        begin_ptr = end_ptr + 1;
      }

      if (num_values != 12) {
        throw std::invalid_argument(
            "[bioparser::MhapParser] error: invalid file format");
      }

      dst.emplace_back(std::unique_ptr<T>(new T(
          lhs_id, rhs_id,
          error,
          num_minmers,
          lhs_strand, lhs_begin, lhs_end, lhs_len,
          rhs_strand, rhs_begin, rhs_end, rhs_len)));

      parsed_bytes += this->storage_ptr();
      this->Clear();
    };

    bool is_eof = false;

    while (true) {
      auto buffer_ptr = this->buffer_ptr();
      for (; buffer_ptr < this->buffer_bytes(); ++buffer_ptr) {
        auto c = this->buffer()[buffer_ptr];
        if (c == '\n') {
          this->Store(buffer_ptr - this->buffer_ptr());
          create_T();
          if (parsed_bytes >= bytes) {
            return dst;
          }
        }
      }
      if (this->buffer_ptr() < buffer_ptr) {
        this->Store(buffer_ptr - this->buffer_ptr());
      }

      if (is_eof) {
        break;
      }
      is_eof = this->Read();
    }

    if (this->storage_ptr() != 0) {
      create_T();
    }

    return dst;
  }

 private:
  explicit MhapParser(gzFile file)
      : Parser<T>(file, 65536) {}  // 64 kB

  friend Parser<T>;
};

}  // namespace bioparser

#endif  // BIOPARSER_MHAP_PARSER_HPP_
