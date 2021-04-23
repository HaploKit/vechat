// Copyright (c) 2020 Robert Vaser

#include "bioparser/fastq_parser.hpp"

#include <numeric>

#include "biosoup/sequence.hpp"
#include "gtest/gtest.h"

namespace bioparser {
namespace test {

class BioparserFastqTest: public ::testing::Test {
 public:
  void Setup(const std::string& file) {
    p = Parser<biosoup::Sequence>::Create<FastqParser>(BIOPARSER_DATA_PATH + file);  // NOLINT
  }

  void Check() {
    EXPECT_EQ(13, s.size());
    EXPECT_EQ(17, std::accumulate(s.begin(), s.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<biosoup::Sequence>& it) {
          return s + it->name.size();
        }));
    EXPECT_EQ(108140, std::accumulate(s.begin(), s.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<biosoup::Sequence>& it) {
          return s + it->data.size();
        }));
    EXPECT_EQ(108140, std::accumulate(s.begin(), s.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<biosoup::Sequence>& it) {
          return s + it->quality.size();
        }));
  }

  std::unique_ptr<Parser<biosoup::Sequence>> p;
  std::vector<std::unique_ptr<biosoup::Sequence>> s;
};

TEST_F(BioparserFastqTest, ParseWhole) {
  Setup("sample.fastq");
  s = p->Parse(-1);
  Check();
}

TEST_F(BioparserFastqTest, ParseInChunks) {
  Setup("sample.fastq");
  for (auto t = p->Parse(65536); !t.empty(); t = p->Parse(65536)) {
    s.insert(
        s.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserFastqTest, FormatError) {
  Setup("sample.fasta");
  try {
    s = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::FastqParser] error: invalid file format");
  }
}

TEST_F(BioparserFastqTest, CompressedParseWhole) {
  Setup("sample.fastq.gz");
  s = p->Parse(-1);
  Check();
}

TEST_F(BioparserFastqTest, CompressedParseInChunks) {
  Setup("sample.fastq.gz");
  for (auto t = p->Parse(65536); !t.empty(); t = p->Parse(65536)) {
    s.insert(
        s.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserFastqTest, CompressedFormatError) {
  Setup("sample.fasta.gz");
  try {
    s = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::FastqParser] error: invalid file format");
  }
}

}  // namespace test
}  // namespace bioparser
