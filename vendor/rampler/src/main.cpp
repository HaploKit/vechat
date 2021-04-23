// Copyright (c) 2021 Robert Vaser

#include <getopt.h>

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

#include "sampler.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace {

static const char* rampler_version = RAMPLER_VERSION;

static struct option options[] = {
  {"out-directory", required_argument, nullptr, 'o'},
  {"version", no_argument, nullptr, 'v'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path, std::string* name, std::string* ext) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  std::size_t c = path.rfind('/');
  *name = (c == std::string::npos ? path : path.substr(c + 1));

  c = name->find('.');
  *name = (c == std::string::npos ? *name : name->substr(0, c));

  if (is_suffix(path, ".fasta")    || is_suffix(path, ".fa") ||
      is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
    *ext = ".fasta";
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq")    || is_suffix(path, ".fq") ||
      is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
    *ext = ".fastq";
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[rampler::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: rampler [options ...] <mode>\n"
      "\n"
      "  <mode>\n"
      "    subsample <sequences> <reference length> <coverage> [<coverage> ...]\n"  // NOLINT
      "\n"
      "      <sequences>\n"
      "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "      <reference length>\n"
      "        integer denoting length of the reference genome (or assembly)\n"
      "      <coverage>\n"
      "        integer denoting desired coverage of the subsampled sequences\n"
      "\n"
      "    split <sequences> <chunk size>\n"
      "\n"
      "      <sequences>\n"
      "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "      <chunk size>\n"
      "        integer denoting the desired chunk size in bytes\n"
      "\n"
      "  options:\n"
      "    -o, --out-directory <string>\n"
      "      default: current directory\n"
      "      path in which sampled files will be created\n"
      "    --version\n"
      "      prints the version number\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> input_parameters;
  std::string out_directory = ".";

  int arg;
  while ((arg = getopt_long(argc, argv, "o:h", options, nullptr)) != -1) {
    switch (arg) {
      case 'o': out_directory = optarg; break;
      case 'v': std::cerr << rampler_version << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  if (optind == argc) {
    std::cerr << "[rampler::] error: missing arguments!" << std::endl;
    return 1;
  }

  for (std::int32_t i = optind; i < argc; ++i) {
    input_parameters.emplace_back(argv[i]);
  }

  bool subsample = false, split = false;
  if (input_parameters[0] == "subsample") {
    subsample = true;
  } else if (input_parameters[0] == "split") {
    split = true;
  } else {
    std::cerr << "[rampler::] error: unknown mode!" << std::endl;
    return 1;
  }

  if ((subsample && input_parameters.size() < 4) ||
      (split && input_parameters.size() < 3)) {
    std::cerr << "[rampler::] error: missing arguments!" << std::endl;
    return 1;
  }

  std::string name{}, ext{};
  auto sparser = CreateParser(input_parameters[1], &name, &ext);
  if (sparser == nullptr) {
    return 1;
  }

  rampler::Sampler sampler{std::move(sparser), name, ext};
  sampler.Initialize();

  if (split) {
    sampler.Split(out_directory, atoi(input_parameters[2].c_str()));
  } else if (subsample) {
    std::uint32_t reference_length = atoi(input_parameters[2].c_str());
    for (std::uint32_t i = 3; i < input_parameters.size(); ++i) {
      sampler.Subsample(
          out_directory,
          reference_length,
          atoi(input_parameters[i].c_str()));
    }
  }

  return 0;
}
