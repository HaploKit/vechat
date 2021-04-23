// Copyright (c) 2020 Robert Vaser

#include "bioparser/sam_parser.hpp"

#include <numeric>
#include <string>

#include "biosoup/overlap.hpp"
#include "gtest/gtest.h"

namespace bioparser {
namespace test {

struct SamOverlap: public biosoup::Overlap {
 public:
  SamOverlap(
      const char* q_name, std::uint32_t q_name_len,
      std::uint32_t flag,
      const char* t_name, std::uint32_t t_name_len,
      std::uint32_t t_begin,
      std::uint32_t map_quality,
      const char* cigar, std::uint32_t cigar_len,
      const char* t_next_name, std::uint32_t t_next_name_len,
      std::uint32_t t_next_begin,
      std::uint32_t template_len,
      const char* data, std::uint32_t data_len,
      const char* quality, std::uint32_t quality_len)
      : biosoup::Overlap(
          0, 0, 0,
          0, t_begin, 0,
          0,
          cigar, cigar_len),
        q_name(q_name, q_name_len),
        flag(flag),
        t_name(t_name, t_name_len),
        map_quality(map_quality),
        t_next_name(t_next_name, t_next_name_len),
        t_next_begin(t_next_begin),
        template_len(template_len),
        data(data, data_len),
        quality(quality, quality_len) {}

  std::string q_name;
  std::uint32_t flag;
  std::string t_name;
  std::uint32_t map_quality;
  std::string t_next_name;
  std::uint32_t t_next_begin;
  std::uint32_t template_len;
  std::string data;
  std::string quality;
};

class BioparserSamTest: public ::testing::Test {
 public:
  void Setup(const std::string& file) {
    p = Parser<SamOverlap>::Create<SamParser>(BIOPARSER_DATA_PATH + file);
  }

  void Check() {
    EXPECT_EQ(48, o.size());
    EXPECT_EQ(795237, std::accumulate(o.begin(), o.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<SamOverlap>& it) {
          return s +
              it->q_name.size() +
              it->t_name.size() +
              it->alignment.size() +
              it->t_next_name.size() +
              it->data.size() +
              it->quality.size();
        }));
    EXPECT_EQ(639677, std::accumulate(o.begin(), o.end(), 0,
        [] (std::uint32_t s, const std::unique_ptr<SamOverlap>& it) {
          return s +
              it->flag +
              it->rhs_begin +
              it->map_quality +
              it->t_next_begin +
              it->template_len;
        }));
  }

  std::unique_ptr<Parser<SamOverlap>> p;
  std::vector<std::unique_ptr<SamOverlap>> o;
};

TEST_F(BioparserSamTest, ParseWhole) {
  Setup("sample.sam");
  o = p->Parse(-1);
  Check();
}

TEST_F(BioparserSamTest, ParseInChunks) {
  Setup("sample.sam");
  for (auto t = p->Parse(1024); !t.empty(); t = p->Parse(1024)) {
    o.insert(
        o.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserSamTest, FormatError) {
  Setup("sample.paf");
  try {
    o = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::SamParser] error: invalid file format");
  }
}

TEST_F(BioparserSamTest, CompressedParseWhole) {
  Setup("sample.sam.gz");
  o = p->Parse(-1);
  Check();
}

TEST_F(BioparserSamTest, CompressedParseInChunks) {
  Setup("sample.sam.gz");
  for (auto t = p->Parse(1024); !t.empty(); t = p->Parse(1024)) {
    o.insert(
        o.end(),
        std::make_move_iterator(t.begin()),
        std::make_move_iterator(t.end()));
  }
  Check();
}

TEST_F(BioparserSamTest, CompressedFormatError) {
  Setup("sample.paf.gz");
  try {
    o = p->Parse(-1);
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[bioparser::SamParser] error: invalid file format");
  }
}

}  // namespace test
}  // namespace bioparser
