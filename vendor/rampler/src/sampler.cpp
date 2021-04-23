// Copyright (c) 2021 Robert Vaser

#include "sampler.hpp"

#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <utility>

namespace rampler {

constexpr uint32_t kChunkSize = 1024 * 1024 * 1024;  // ~ 1GB

Sampler::Sampler(
    std::unique_ptr<bioparser::Parser<biosoup::Sequence>> sparser,
    const std::string& base_name,
    const std::string& extension)
    : sparser_(std::move(sparser)),
      sequences_length_(0),
      base_name_(base_name),
      extension_(extension) {
}

void Sampler::Initialize() {
  if (sequences_length_ != 0) {
    return;
  }

  sparser_->Reset();
  while (true) {
    auto sequences = sparser_->Parse(1ULL << 30);
    if (sequences.empty()) {
      break;
    }

    for (const auto& it : sequences) {
      sequences_length_ += it->data.size();
    }
  }
}

void Sampler::Subsample(
    const std::string& out_directory,
    std::uint32_t reference_length,
    std::uint32_t coverage) {
  if (coverage * reference_length > sequences_length_) {
    std::cerr << "[rampler::Sampler::subsample] warning: "
              << "insufficient data for coverage of " << coverage
              << std::endl;
    return;
  }

  std::random_device r;
  std::mt19937 generator(r());
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  double ratio = (coverage * reference_length) / static_cast<double>(sequences_length_);  // NOLINT

  std::string out_path =
      out_directory + "/" + base_name_ + "_" +
      std::to_string(coverage) + "x" + extension_;

  std::ofstream ofs(out_path);
  if (!ofs.is_open()) {
    throw std::runtime_error(
        "[rampler::Sampler::subsample] error: unable to create file on disk!");
  }

  sparser_->Reset();
  while (true) {
    auto sequences = sparser_->Parse(1ULL << 30);
    if (sequences.empty()) {
      break;
    }

    for (const auto& it : sequences) {
      if (distribution(generator) < ratio) {
        if (it->quality.empty()) {
          ofs << ">" << it->name << std::endl
              << it->data << std::endl;
        } else {
          ofs << "@" << it->name << std::endl
              << it->data << std::endl
              << "+" << std::endl
              << it->quality << std::endl;
        }
      }
    }
  }

  ofs.close();
}

void Sampler::Split(const std::string& out_directory, std::uint32_t chunk_size) {  // NOLINT
  if (chunk_size > sequences_length_) {
    std::cerr << "[rampler::Sampler::split] warning: "
              << "insufficient data for chunk size " << chunk_size
              << std::endl;
    return;
  }

  uint32_t chunk_number = 0;

  sparser_->Reset();
  while (true) {
    auto sequences = sparser_->Parse(chunk_size);
    if (sequences.empty()) {
      break;
    }

    std::string out_path =
        out_directory + "/" + base_name_ + "_" +
        std::to_string(chunk_number++) + extension_;

    std::ofstream ofs(out_path);
    if (!ofs.is_open()) {
      throw std::runtime_error(
          "[rampler::Sampler::subsample] error: unable to create file on disk!");  // NOLINT
    }

    for (const auto& it : sequences) {
      if (it->quality.empty()) {
        ofs << ">" << it->name << std::endl
            << it->data << std::endl;
      } else {
        ofs << "@" << it->name << std::endl
            << it->data << std::endl
            << "+" << std::endl
            << it->quality << std::endl;
      }
    }

    ofs.close();
  }
}

}  // namespace rampler
