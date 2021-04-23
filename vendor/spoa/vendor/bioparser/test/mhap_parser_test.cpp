// Copyright (c) 2020 Robert Vaser

#include "bioparser/mhap_parser.hpp"

#include <numeric>
#include <string>

#include "biosoup/overlap.hpp"
#include "gtest/gtest.h"

namespace bioparser {
namespace test {

struct MhapOverlap: public biosoup::Overlap {
 public:
  MhapOverlap(
      std::uint64_t lhs_id,
      std::uint64_t rhs_id,
      double error,
      std::uint32_t num_minmers,
      std::uint32_t lhs_strand,
      std::uint32_t lhs_begin,
      std::uint32_t lhs_end,
      std::uint32_t lhs_len,
      std::uint32_t rhs_strand,
      std::uint32_t rhs_begin,
      std::uint32_t rhs_end,
      std::uint32_t rhs_len)
      : biosoup::Overlap(
          lhs_id, lhs_begin, lhs_end,
          rhs_id, rhs_begin, rhs_end,
          num_minmers,
          lhs_strand == rhs_strand),
        error(error * 10000),
        lhs_len(lhs_len),
        rhs_len(rhs_len) {}

  std::uint32_t error;
  std::uint32_t lhs_len;
  std::uint32_t rhs_len;
};

class BioparserMhapTest: public ::testing::Test {
 public:
  void Setup(const std::string& file) {
    p = Parser<MhapOverlap>::Create<MhapParser>(BIOPARSER_DATA_PATH + file);
  }

  void Check() {
    EXPECT_EQ(150, o.size());
    EXPECT_EQ(7816660, std::accumulate(o.begin(), o.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<MhapOverlap>& it) {
          return s +
              it->lhs_id + it->lhs_begin + it->lhs_end + it->lhs_len +
              it->rhs_id + it->rhs_begin + it->rhs_end + it->rhs_len +
              it->score +
              it->strand +
              it->error;
        }));
  }

  std::unique_ptr<Parser<MhapOverlap>> p;
  std::vector<std::unique_ptr<MhapOverlap>> o;
};

TEST_F(BioparserMhapTest, ParseWhole) {
  Setup("sample.mhap");
  o = p->Parse(-1);
  Check();
}

TEST_F(BioparserMhapTest, ParseInChunks) {
  Setup("sample.mhap");
  for (auto t = p->Parse(1024); !t.empty(); t = p->Parse(1024)) {
    o.insert(
        o.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserMhapTest, FormatError) {
  Setup("sample.paf");
  try {
    o = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::MhapParser] error: invalid file format");
  }
}

TEST_F(BioparserMhapTest, CompressedParseWhole) {
  Setup("sample.mhap.gz");
  o = p->Parse(-1);
  Check();
}

TEST_F(BioparserMhapTest, CompressedParseInChunks) {
  Setup("sample.mhap.gz");
  for (auto t = p->Parse(1024); !t.empty(); t = p->Parse(1024)) {
    o.insert(
        o.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserMhapTest, CompressedFormatError) {
  Setup("sample.paf.gz");
  try {
    o = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::MhapParser] error: invalid file format");
  }
}

}  // namespace test
}  // namespace bioparser
