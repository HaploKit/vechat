// Copyright (c) 2020 Robert Vaser

#ifndef BIOPARSER_FASTQ_PARSER_HPP_
#define BIOPARSER_FASTQ_PARSER_HPP_

#include <cstdint>
#include <memory>
#include <vector>
#include <stdexcept>
#include <string>

#include "bioparser/parser.hpp"

namespace bioparser {

template<class T>
class FastqParser: public Parser<T> {
 public:
  FastqParser(const FastqParser&) = delete;
  FastqParser& operator=(const FastqParser&) = delete;

  FastqParser(FastqParser&&) = delete;
  FastqParser& operator=(FastqParser&&) = delete;

  ~FastqParser() {}

  std::vector<std::unique_ptr<T>> Parse(
      std::uint64_t bytes, bool shorten_names = true) override {
    std::vector<std::unique_ptr<T>> dst;
    std::uint64_t parsed_bytes = 0;
    std::uint32_t data_ptr = 0;
    std::uint32_t comment_ptr = 0;
    std::uint32_t quality_ptr = 0;

    auto create_T = [&] () -> void {
      if (data_ptr == 0 || comment_ptr == 0 || quality_ptr == 0) {
        throw std::invalid_argument(
            "[bioparser::FastqParser] error: invalid file format");
      }

      auto name_len = shorten_names ?
          this->Shorten(this->storage().data(), data_ptr) :
          this->RightStrip(this->storage().data(), data_ptr);

      auto data_len = comment_ptr - data_ptr;

      auto quality_len = this->storage_ptr() - quality_ptr;

      if (name_len == 0 || this->storage()[0] != '@' || data_len == 0 ||
          quality_len == 0 || data_len != quality_len) {
        throw std::invalid_argument(
            "[bioparser::FastqParser] error: invalid file format");
      }

      dst.emplace_back(std::unique_ptr<T>(new T(
          static_cast<const char*>(this->storage().data() + 1), name_len - 1,
          static_cast<const char*>(this->storage().data() + data_ptr), data_len,
          static_cast<const char*>(this->storage().data() + quality_ptr), quality_len)));  // NOLINT

      parsed_bytes += this->storage_ptr();
      data_ptr = 0;
      comment_ptr = 0;
      quality_ptr = 0;
      this->Clear();
    };

    bool is_eof = false;
    bool is_name = true;
    bool is_data = false;
    bool is_comment = false;
    bool is_quality = false;

    while (true) {
      auto buffer_ptr = this->buffer_ptr();
      for (; buffer_ptr < this->buffer_bytes(); ++buffer_ptr) {
        auto c = this->buffer()[buffer_ptr];
        if (c == '\n') {
          this->Store(buffer_ptr - this->buffer_ptr(), !is_name);
          if (is_name) {
            is_name = false;
            is_data = true;
            data_ptr = this->storage_ptr();
          } else if (is_comment) {
            is_comment = false;
            is_quality = true;
            quality_ptr = this->storage_ptr();
          } else if (is_quality &&
                this->storage_ptr() - quality_ptr == comment_ptr - data_ptr) {
            is_quality = false;
            is_name = true;
            create_T();
            if (parsed_bytes >= bytes) {
              return dst;
            }
          }
        } else if (is_data && c == '+') {
          is_data = false;
          is_comment = true;
          comment_ptr = this->storage_ptr();
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
  explicit FastqParser(gzFile file)
      : Parser<T>(file, 4194304) {}  // 4 MB

  friend Parser<T>;
};

}  // namespace bioparser

#endif  // BIOPARSER_FASTQ_PARSER_HPP_
