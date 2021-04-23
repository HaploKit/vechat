// Copyright (c) 2020 Robert Vaser

#include "bioparser/fasta_parser.hpp"

#include <numeric>

#include "biosoup/sequence.hpp"
#include "gtest/gtest.h"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace bioparser {
namespace test {

class BioparserFastaTest: public ::testing::Test {
 public:
  void Setup(const std::string& file) {
    p = Parser<biosoup::Sequence>::Create<FastaParser>(BIOPARSER_DATA_PATH + file);  // NOLINT
  }

  void Check(bool is_trimmed = true) {
    EXPECT_EQ(14, s.size());
    EXPECT_EQ(65 + !is_trimmed * 10, std::accumulate(s.begin(), s.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<biosoup::Sequence>& it) {
          return s + it->name.size();
        }));
    EXPECT_EQ(109117, std::accumulate(s.begin(), s.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<biosoup::Sequence>& it) {
          return s + it->data.size();
        }));
    EXPECT_EQ(0, std::accumulate(s.begin(), s.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<biosoup::Sequence>& it) {
          return s + it->quality.size();
        }));
  }

  std::unique_ptr<Parser<biosoup::Sequence>> p;
  std::vector<std::unique_ptr<biosoup::Sequence>> s;
};

TEST_F(BioparserFastaTest, ParseWhole) {
  Setup("sample.fasta");
  s = p->Parse(-1);
  Check();
}

TEST_F(BioparserFastaTest, ParseWholeWithoutTrimming) {
  Setup("sample.fasta");
  s = p->Parse(-1, false);
  Check(false);
}

TEST_F(BioparserFastaTest, ParseInChunks) {
  Setup("sample.fasta");
  for (auto t = p->Parse(65536); !t.empty(); t = p->Parse(65536)) {
    s.insert(
        s.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserFastaTest, FormatError) {
  Setup("sample.fastq");
  try {
    s = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::FastaParser] error: invalid file format");
  }
}

TEST_F(BioparserFastaTest, ParseAndReset) {
  Setup("sample.fasta");
  s = p->Parse(-1);
  p->Reset();
  s.clear();
  for (auto t = p->Parse(65536); !t.empty(); t = p->Parse(65536)) {
    s.insert(
        s.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserFastaTest, CompressedParseWhole) {
  Setup("sample.fasta.gz");
  s = p->Parse(-1);
  Check();
}

TEST_F(BioparserFastaTest, CompressedParseWholeWithoutTrimming) {
  Setup("sample.fasta.gz");
  s = p->Parse(-1, false);
  Check(false);
}

TEST_F(BioparserFastaTest, CompressedParseInChunks) {
  Setup("sample.fasta.gz");
  for (auto t = p->Parse(65536); !t.empty(); t = p->Parse(65536)) {
    s.insert(
        s.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserFastaTest, CompressedParseAndReset) {
  Setup("sample.fasta.gz");
  s = p->Parse(-1);
  p->Reset();
  s.clear();
  for (auto t = p->Parse(65536); !t.empty(); t = p->Parse(65536)) {
    s.insert(
        s.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserFastaTest, CompressedFormatError) {
  Setup("sample.fastq.gz");
  try {
    s = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::FastaParser] error: invalid file format");
  }
}

}  // namespace test
}  // namespace bioparser
