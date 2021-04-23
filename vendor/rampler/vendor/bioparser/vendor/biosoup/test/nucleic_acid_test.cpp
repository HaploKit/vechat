// Copyright (c) 2020 Robert Vaser

#include "biosoup/nucleic_acid.hpp"

#include <iostream>

#include "gtest/gtest.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace biosoup {
namespace test {

TEST(BiosoupNucleicAcidTest, NucleotideError) {
  try {
    NucleicAcid s{"test", "EFIJLOPQUXZ"};
  } catch (std::invalid_argument& exception) {
    EXPECT_STREQ(
        exception.what(),
        "[biosoup::NucleicAcid::NucleicAcid] error: not a nucleotide");
  }
}

TEST(BiosoupNucleicAcidTest, Inflate) {
  NucleicAcid s{"test", "AaAaCcCcGgGgTtTt------ACGTRYKMSWBDHVN-nvhdbwsmkyrtgca------tTtTgGgGcCcCaAaA"};  // NOLINT
  EXPECT_EQ(s.deflated_data.size(), s.deflated_data.capacity());
  EXPECT_EQ(0, s.Code(16));
  EXPECT_EQ(1, s.Code(23));
  EXPECT_EQ(2, s.Code(35));
  EXPECT_EQ(3, s.Code(59));
  EXPECT_EQ("AAAACCCCGGGGTTTTAAAAAAACGTATGCCACATGAAAGTACACCGTATGCAAAAAAATTTTGGGGCCCCAAAA", s.InflateData());  // NOLINT
  EXPECT_EQ("TATGCCACATGAAAGTACACCGTAT", s.InflateData(25, 25));
  EXPECT_EQ("TGAAAGT", s.InflateData(34, 7));
  EXPECT_EQ("", s.InflateData(75, 42));
  EXPECT_EQ("C", s.InflateData(29, 1));
  EXPECT_EQ("G", s.InflateData(64, 1));
  EXPECT_EQ("CCAAAA", s.InflateData(69));
}

TEST(BiosoupNucleicAcidTest, Quality) {
  NucleicAcid s{
      "test",
      "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",  // NOLINT
      "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"};  // NOLINT
  EXPECT_EQ(s.block_quality.size(), s.block_quality.capacity());
  EXPECT_EQ(31, s.Score(42));
  EXPECT_EQ(78, s.Score(84));
  EXPECT_EQ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@oooooooooooooooooooooooooooooo", s.InflateQuality());  // NOLINT
  EXPECT_EQ("@@@@@@@@@@@@@@@@", s.InflateQuality(0, 16));
  EXPECT_EQ("@@oo", s.InflateQuality(62, 4));
  EXPECT_EQ("@", s.InflateQuality(63, 1));
  EXPECT_EQ("o", s.InflateQuality(64, 1));
  EXPECT_EQ("o", s.InflateQuality(93));
  EXPECT_EQ("ooo", s.InflateQuality(91, 20));
  EXPECT_EQ("", s.InflateQuality(95));
}

TEST(BiosoupNucleicAcidTest, ReverseAndComplement) {
  NucleicAcid s{
      "test",
      "ACGTACTGAGCTAGTCATCGATGCCAGTCATGCGATCGTACTAGCTGAGACTGATCGCATGCTAGTACGTCA",  // NOLINT
      "0123456789012345678901234567890123456789012345678901234567890123ZZZZZZZZ"};  // NOLINT
  NucleicAcid c{s};
  c.ReverseAndComplement();
  EXPECT_EQ("TGACGTACTAGCATGCGATCAGTCTCAGCTAGTACGATCGCATGACTGGCATCGATGACTAGCTCAGTACGT", c.InflateData());  // NOLINT
  EXPECT_EQ("ZZZZZZZZ4444444444444444444444444444444444444444444444444444444444444444", c.InflateQuality());  // NOLINT
  c.ReverseAndComplement();
  EXPECT_EQ(c.InflateData(), s.InflateData());
  EXPECT_EQ(c.InflateQuality(), s.InflateQuality());
}

}  // namespace test
}  // namespace biosoup
