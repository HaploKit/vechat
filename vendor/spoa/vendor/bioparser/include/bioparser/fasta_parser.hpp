// Copyright (c) 2020 Robert Vaser

#ifndef BIOPARSER_FASTA_PARSER_HPP_
#define BIOPARSER_FASTA_PARSER_HPP_

#include <cstdint>
#include <memory>
#include <vector>
#include <stdexcept>
#include <string>

#include "bioparser/parser.hpp"

namespace bioparser {

template<class T>
class FastaParser: public Parser<T> {
 public:
  FastaParser(const FastaParser&) = delete;
  FastaParser& operator=(const FastaParser&) = delete;

  FastaParser(FastaParser&&) = delete;
  FastaParser& operator=(FastaParser&&) = delete;

  ~FastaParser() {}

  std::vector<std::unique_ptr<T>> Parse(
      std::uint64_t bytes, bool shorten_names = true) override {
    std::vector<std::unique_ptr<T>> dst;
    std::uint64_t parsed_bytes = 0;
    std::uint32_t data_ptr = 0;

    auto create_T = [&] () -> void {
      if (data_ptr == 0) {
        throw std::invalid_argument(
            "[bioparser::FastaParser] error: invalid file format");
      }

      auto name_len = shorten_names ?
          this->Shorten(this->storage().data(), data_ptr) :
          this->RightStrip(this->storage().data(), data_ptr);

      auto data_len = this->storage_ptr() - data_ptr;

      if (name_len == 0 || this->storage()[0] != '>' || data_len == 0) {
        throw std::invalid_argument(
            "[bioparser::FastaParser] error: invalid file format");
      }

      dst.emplace_back(std::unique_ptr<T>(new T(
          static_cast<const char*>(this->storage().data() + 1), name_len - 1,
          static_cast<const char*>(this->storage().data() + data_ptr), data_len)));  // NOLINT

      parsed_bytes += this->storage_ptr();
      data_ptr = 0;
      this->Clear();
    };

    bool is_eof = false;
    bool is_name = true;

    while (true) {
      auto buffer_ptr = this->buffer_ptr();
      for (; buffer_ptr < this->buffer_bytes(); ++buffer_ptr) {
        auto c = this->buffer()[buffer_ptr];
        if (c == '\n') {
          this->Store(buffer_ptr - this->buffer_ptr(), !is_name);
          if (is_name) {
            data_ptr = this->storage_ptr();
            is_name = false;
          }
        } else if (!is_name && c == '>') {
          is_name = true;
          create_T();
          if (parsed_bytes >= bytes) {
            return dst;
          }
        }
      }
      if (this->buffer_ptr() < buffer_ptr) {
        this->Store(buffer_ptr - this->buffer_ptr(), !is_name);
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
  explicit FastaParser(gzFile file)
      : Parser<T>(file, 4194304) {}  // 4 MB

  friend Parser<T>;
};

}  // namespace bioparser

#endif  // BIOPARSER_FASTA_PARSER_HPP_
