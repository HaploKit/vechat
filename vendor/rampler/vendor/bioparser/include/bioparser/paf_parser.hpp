// Copyright (c) 2020 Robert Vaser

#ifndef BIOPARSER_PAF_PARSER_HPP_
#define BIOPARSER_PAF_PARSER_HPP_

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <vector>
#include <stdexcept>

#include "bioparser/parser.hpp"

namespace bioparser {

template<class T>
class PafParser: public Parser<T> {
 public:
  PafParser(const PafParser&) = delete;
  PafParser& operator=(const PafParser&) = delete;

  PafParser(PafParser&&) = delete;
  PafParser& operator=(PafParser&&) = delete;

  ~PafParser() {}

  std::vector<std::unique_ptr<T>> Parse(
      std::uint64_t bytes, bool shorten_names = true) override {
    std::vector<std::unique_ptr<T>> dst;
    std::uint64_t parsed_bytes = 0;

    const char* q_name = nullptr;
    std::uint32_t q_name_len = 0;
    std::uint32_t q_len = 0;
    std::uint32_t q_begin = 0;
    std::uint32_t q_end = 0;
    const char* t_name = nullptr;
    std::uint32_t t_name_len = 0;
    std::uint32_t t_len = 0;
    std::uint32_t t_begin = 0;
    std::uint32_t t_end = 0;
    std::uint32_t num_matches = 0;
    std::uint32_t overlap_len = 0;
    std::uint32_t quality = 0;
    char orientation = '\0';

    auto create_T = [&] () -> void {
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
          case 1: q_len = std::atoi(this->storage().data() + begin_ptr); break;
          case 2: q_begin = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 3: q_end = std::atoi(this->storage().data() + begin_ptr); break;
          case 4: orientation = this->storage()[begin_ptr]; break;
          case 5:
            t_name = this->storage().data() + begin_ptr;
            t_name_len = end_ptr - begin_ptr;
            break;
          case 6: t_len = std::atoi(this->storage().data() + begin_ptr); break;
          case 7: t_begin = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 8: t_end = std::atoi(this->storage().data() + begin_ptr); break;
          case 9: num_matches = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 10: overlap_len = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
          case 11: quality = std::atoi(this->storage().data() + begin_ptr); break;  // NOLINT
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
            "[bioparser::PafParser] error: invalid file format");
      }

      q_name_len = shorten_names ?
          this->Shorten(q_name, q_name_len) :
          this->RightStrip(q_name, q_name_len);

      t_name_len = shorten_names ?
          this->Shorten(t_name, t_name_len) :
          this->RightStrip(t_name, t_name_len);

      if (q_name_len == 0 || t_name_len == 0) {
        throw std::invalid_argument(
            "[bioparser::PafParser] error: invalid file format");
      }

      dst.emplace_back(std::unique_ptr<T>(new T(
          q_name, q_name_len, q_len, q_begin, q_end,
          orientation,
          t_name, t_name_len, t_len, t_begin, t_end,
          num_matches,
          overlap_len,
          quality)));

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
  explicit PafParser(gzFile file)
      : Parser<T>(file, 65536) {}  // 64 kB

  friend Parser<T>;
};

}  // namespace bioparser

#endif  // BIOPARSER_PAF_PARSER_HPP_
