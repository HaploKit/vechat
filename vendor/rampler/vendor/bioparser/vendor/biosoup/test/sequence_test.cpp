// Copyright (c) 2020 Robert Vaser

#include "biosoup/sequence.hpp"

#include "gtest/gtest.h"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace biosoup {
namespace test {

TEST(BiosoupSequenceTest, ReverseAndComplement) {
  Sequence s{"Test", "ACGTURYKMSWBDHVN", "0123456789:;<=>?"};
  s.ReverseAndComplement();
  EXPECT_EQ(0, s.id);
  EXPECT_EQ("NBDHVWSKMRYAACGT", s.data);
  EXPECT_EQ("?>=<;:9876543210", s.quality);
}

}  // namespace test
}  // namespace biosoup
