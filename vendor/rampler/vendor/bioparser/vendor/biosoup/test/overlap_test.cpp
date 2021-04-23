// Copyright (c) 2020 Robert Vaser

#include "biosoup/overlap.hpp"

#include "gtest/gtest.h"

namespace biosoup {
namespace test {

TEST(BiosoupOverlapTest, Compile) {
  Overlap o{0, 10, 20, 1, 0, 10, 10, 0};
}

}  // namespace test
}  // namespace biosoup
