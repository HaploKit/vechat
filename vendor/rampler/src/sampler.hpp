// Copyright (c) 2021 Robert Vaser

#ifndef RAMPLER_SAMPLER_HPP_
#define RAMPLER_SAMPLER_HPP_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "bioparser/parser.hpp"
#include "biosoup/sequence.hpp"

namespace rampler {

class Sampler;
std::unique_ptr<Sampler> createSampler(const std::string& sequences_path);

class Sampler {
 public:
  Sampler(
      std::unique_ptr<bioparser::Parser<biosoup::Sequence>> sparser,
      const std::string& base_name,
      const std::string& extension);

  Sampler(const Sampler&) = delete;
  Sampler& operator=(const Sampler&) = delete;

  Sampler(Sampler&&) = delete;
  Sampler& operator=(Sampler&&) = delete;

  ~Sampler() = default;

  void Initialize();

  void Subsample(
      const std::string& out_directory,
      std::uint32_t reference_length,
      std::uint32_t coverage);

  void Split(const std::string& out_directory, std::uint32_t chunk_size);

 private:
  std::unique_ptr<bioparser::Parser<biosoup::Sequence>> sparser_;
  std::uint64_t sequences_length_;
  std::string base_name_;
  std::string extension_;
};

}  // namespace rampler

#endif  // RAMPLER_SAMPLER_HPP_
