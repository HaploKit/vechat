// Copyright (c) 2020 Robert Vaser

#ifndef BIOPARSER_PARSER_HPP_
#define BIOPARSER_PARSER_HPP_

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "zlib.h"  // NOLINT

namespace bioparser {

template<class T>
class Parser {  // Parser factory
 public:
  Parser(const Parser&) = delete;
  Parser& operator=(const Parser&) = delete;

  Parser(Parser&&) = delete;
  Parser& operator=(Parser&&) = delete;

  virtual ~Parser() {}

  template<template<class> class P>
  static std::unique_ptr<Parser<T>> Create(const std::string& path) {
    auto file = gzopen(path.c_str(), "r");
    if (file == nullptr) {
      throw std::invalid_argument(
          "[bioparser::Parser::Create] error: unable to open file " + path);
    }
    return std::unique_ptr<Parser<T>>(new P<T>(file));
  }

  // by default, all parsers shrink sequence names to the first white space
  virtual std::vector<std::unique_ptr<T>> Parse(
      std::uint64_t bytes, bool shorten_names = true) = 0;

  void Reset() {
    gzseek(file_.get(), 0, SEEK_SET);
    buffer_ptr_ = 0;
    buffer_bytes_ = 0;
  }

 protected:
  Parser(gzFile file, std::uint32_t storage_size)
      : file_(file, gzclose),
        buffer_(65536, 0),  // 64 kB
        buffer_ptr_(0),
        buffer_bytes_(0),
        storage_(storage_size, 0),
        storage_ptr_(0) {}

  const std::vector<char>& buffer() const {
    return buffer_;
  }

  std::uint32_t buffer_ptr() const {
    return buffer_ptr_;
  }

  std::uint32_t buffer_bytes() const {
    return buffer_bytes_;
  }

  const std::vector<char>& storage() const {
    return storage_;
  }

  std::uint32_t storage_ptr() const {
    return storage_ptr_;
  }

  bool Read() {
    buffer_ptr_ = 0;
    buffer_bytes_ = gzread(file_.get(), buffer_.data(), buffer_.size());
    return buffer_bytes_ < buffer_.size();
  }

  void Store(std::uint32_t count, bool strip = false) {
    if (buffer_ptr_ + count > buffer_.size()) {
      throw std::invalid_argument(
          "[bioparser::Parser::Store] error: buffer overflow");
    }
    if (storage_ptr_ + count > storage_.size()) {
      storage_.resize(2 * storage_.size());
    }
    std::memcpy(&storage_[storage_ptr_], &buffer_[buffer_ptr_], count);
    storage_ptr_ += strip ? RightStrip(&storage_[storage_ptr_], count) : count;
    buffer_ptr_ += count + 1;  // ignore sought character
  }

  void Terminate(std::uint32_t i) {
    storage_[i] = '\0';
  }

  void Clear() {
    storage_ptr_ = 0;
  }

  static std::uint32_t RightStrip(const char* str, std::uint32_t str_len) {
    while (str_len > 0 && std::isspace(str[str_len - 1])) {
      --str_len;
    }
    return str_len;
  }

  static std::uint32_t Shorten(const char* str, std::uint32_t str_len) {
    for (std::uint32_t i = 0; i < str_len; ++i) {
      if (std::isspace(str[i])) {
        return i;
      }
    }
    return str_len;
  }

 private:
  std::unique_ptr<gzFile_s, int(*)(gzFile)> file_;
  std::vector<char> buffer_;
  std::uint32_t buffer_ptr_;
  std::uint32_t buffer_bytes_;
  std::vector<char> storage_;
  std::uint32_t storage_ptr_;
};

}  // namespace bioparser

#endif  // BIOPARSER_PARSER_HPP_
