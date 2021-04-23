// Copyright (c) 2020 Robert Vaser

#include "bioparser/paf_parser.hpp"

#include <numeric>
#include <string>

#include "biosoup/overlap.hpp"
#include "gtest/gtest.h"

namespace bioparser {
namespace test {

struct PafOverlap: public biosoup::Overlap {
 public:
  PafOverlap(
      const char* q_name, std::uint32_t q_name_len,
      std::uint32_t q_len,
      std::uint32_t q_begin,
      std::uint32_t q_end,
      char orientation,
      const char* t_name, std::uint32_t t_name_len,
      std::uint32_t t_len,
      std::uint32_t t_begin,
      std::uint32_t t_end,
      std::uint32_t score,
      std::uint32_t overlap_len,
      std::uint32_t quality)
      : biosoup::Overlap(
          0, q_begin, q_end,
          0, t_begin, t_end,
          score,
          orientation == '+'),
        q_name(q_name, q_name_len),
        q_len(q_len),
        t_name(t_name, t_name_len),
        t_len(t_len),
        overlap_len(overlap_len),
        quality(quality) {}

  std::string q_name;
  std::uint32_t q_len;
  std::string t_name;
  std::uint32_t t_len;
  std::uint32_t overlap_len;
  std::uint32_t quality;
};

class BioparserPafTest: public ::testing::Test {
 public:
  void Setup(const std::string& file) {
    p = Parser<PafOverlap>::Create<PafParser>(BIOPARSER_DATA_PATH + file);
  }

  void Check() {
    EXPECT_EQ(500, o.size());
    EXPECT_EQ(96478, std::accumulate(o.begin(), o.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<PafOverlap>& it) {
          return s + it->q_name.size() + it->t_name.size();
        }));
    EXPECT_EQ(18472506, std::accumulate(o.begin(), o.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<PafOverlap>& it) {
          return s +
              it->lhs_id + it->lhs_begin + it->lhs_end + it->q_len +
              it->rhs_id + it->rhs_begin + it->rhs_end + it->t_len +
              it->score +
              it->strand +
              it->overlap_len +
              it->quality;
        }));
  }

  std::unique_ptr<Parser<PafOverlap>> p;
  std::vector<std::unique_ptr<PafOverlap>> o;
};

TEST_F(BioparserPafTest, ParseWhole) {
  Setup("sample.paf");
  o = p->Parse(-1);
  Check();
}

TEST_F(BioparserPafTest, ParseInChunks) {
  Setup("sample.paf");
  for (auto t = p->Parse(1024); !t.empty(); t = p->Parse(1024)) {
    o.insert(
        o.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserPafTest, FormatError) {
  Setup("sample.mhap");
  try {
    o = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::PafParser] error: invalid file format");
  }
}

TEST_F(BioparserPafTest, CompressedParseWhole) {
  Setup("sample.paf.gz");
  o = p->Parse(-1);
  Check();
}

TEST_F(BioparserPafTest, CompressedParseInChunks) {
  Setup("sample.paf.gz");
  for (auto t = p->Parse(1024); !t.empty(); t = p->Parse(1024)) {
    o.insert(
        o.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserPafTest, CompressedFormatError) {
  Setup("sample.mhap.gz");
  try {
    o = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::PafParser] error: invalid file format");
  }
}

}  // namespace test
}  // namespace bioparser
