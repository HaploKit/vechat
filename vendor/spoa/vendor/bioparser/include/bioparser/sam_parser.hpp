// Copyright (c) 2020 Robert Vaser

#ifndef BIOPARSER_SAM_PARSER_HPP_
#define BIOPARSER_SAM_PARSER_HPP_

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <vector>
#include <stdexcept>

#include "bioparser/parser.hpp"

namespace bioparser {

template<class T>
class SamParser: public Parser<T> {
 public:
  SamParser(const SamParser&) = delete;
  SamParser& operator=(const SamParser&) = delete;

  SamParser(SamParser&&) = delete;
  SamParser& operator=(SamParser&&) = delete;

  ~SamParser() {}

  std::vector<std::unique_ptr<T>> Parse(
      std::uint64_t bytes, bool shorten_names = true) override {
    std::vector<std::unique_ptr<T>> dst;
    std::uint64_t parsed_bytes = 0;

    const char* q_name = nullptr;
    std::uint32_t q_name_len = 0;
    std::uint32_t flag = 0;
    const char* t_name = nullptr;
    std::uint32_t t_name_len = 0;
    std::uint32_t t_begin = 0;
    std::uint32_t map_quality = 0;
    const char* cigar = nullptr;
    std::uint32_t cigar_len = 0;
    const char* t_next_name = nullptr;
    std::uint32_t t_next_name_len = 0;
    std::uint32_t t_next_begin = 0;
    std::uint32_t template_len = 0;
    const char* data = nullptr;
    std::uint32_t data_len = 0;
    const char* quality = nullptr;
    std::uint32_t quality_len = 0;

    auto create_T = [&] () -> void {
      if (this->storage()[0] == '@') {  // file header
        this->Clear();
        return;
      }
      auto storage_ptr = this->RightStrip(
          this->storage().data(),
          this->storage_ptr());
      this->Terminate(storage_ptr);

      std::uint32_t num_values = 0;
      std::uint32_t begin_ptr = 0;
      while (true) {
        auto end_ptr = begin_ptr;
        while (end_ptr < storage_ptr && this->storage()[end_ptr] != '\t') {
          ++end_ptr;
        }
        this->Terminate(end_ptr);

        switch (num_values) {
          case 0:
            q_name = this->storage().data() + begin_ptr;
            q_name_len = end_ptr - begin_ptr;
            break;
          case 1: flag = std::atoi(this->storage().data() + begin_ptr); break;
          case 2:
            t_name = this->storage().data() + begin_ptr;
            t_name_len = end_ptr - begin_ptr;
            break;
          case 3: t_begin = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 4: map_quality = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 5:
            cigar = this->storage().data() + begin_ptr;
            cigar_len = end_ptr - begin_ptr;
            break;
          case 6:
            t_next_name = this->storage().data() + begin_ptr;
            t_next_name_len = end_ptr - begin_ptr;
            break;
          case 7: t_next_begin = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 8: template_len = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 9:
            data = this->storage().data() + begin_ptr;
            data_len = end_ptr - begin_ptr;
            break;
          case 10:
            quality = this->storage().data() + begin_ptr;
            quality_len = end_ptr - begin_ptr;
            break;
          default: break;
        }

        ++num_values;
        if (end_ptr == storage_ptr || num_values == 11) {
          break;
        }
        begin_ptr = end_ptr + 1;
      }

      if (num_values != 11) {
        throw std::invalid_argument(
            "[bioparser::SamParser] error: invalid file format");
      }

      q_name_len = shorten_names ?
          this->Shorten(q_name, q_name_len) :
          this->RightStrip(q_name, q_name_len);

      t_name_len = shorten_names ?
          this->Shorten(t_name, t_name_len) :
          this->RightStrip(t_name, t_name_len);

      cigar_len = this->RightStrip(cigar, cigar_len);

      t_next_name_len = shorten_names ?
          this->Shorten(t_next_name, t_next_name_len) :
          this->RightStrip(t_next_name, t_next_name_len);

      data_len = this->RightStrip(data, data_len);
      quality_len = this->RightStrip(quality, quality_len);

      if (q_name_len == 0 || t_name_len == 0 || cigar_len == 0 ||
          t_next_name_len == 0 || data_len == 0 || quality_len == 0 ||
          (data_len > 1 && quality_len > 1 && data_len != quality_len)) {
        throw std::invalid_argument(
            "[bioparser::SamParser] error: invalid file format");
      }

      dst.emplace_back(std::unique_ptr<T>(new T(
          q_name, q_name_len,
          flag,
          t_name, t_name_len, t_begin,
          map_quality,
          cigar, cigar_len,
          t_next_name, t_next_name_len, t_next_begin,
          template_len,
          data, data_len,
          quality, quality_len)));

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
  explicit SamParser(gzFile file)
      : Parser<T>(file, 65536) {}  // 64 kB

  friend Parser<T>;
};

}  // namespace bioparser

#endif  // BIOPARSER_SAM_PARSER_HPP_
