// Copyright (c) 2020 Robert Vaser

#include "biosoup/progress_bar.hpp"

#include "gtest/gtest.h"

namespace biosoup {
namespace test {

TEST(BiosoupProgressBarTest, Operator) {
  ProgressBar pb{8, 8};
  EXPECT_EQ("        ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("=>      ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("==>     ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("===>    ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("====>   ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("=====>  ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("======> ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("=======>", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("========", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("========", ::testing::PrintToString(pb));

  pb = {8, 9};
  EXPECT_EQ("        ", ::testing::PrintToString(pb));

  pb = {8, 7};
  EXPECT_EQ("       ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("       ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("=>     ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("==>    ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("===>   ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("====>  ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("=====> ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("======>", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("=======", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("=======", ::testing::PrintToString(pb));

  pb = {8, 4};
  EXPECT_EQ("    ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("    ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("=>  ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("=>  ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("==> ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("==> ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("===>", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("===>", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("====", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("====", ::testing::PrintToString(pb));

  pb = {8, 2};
  EXPECT_EQ("  ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("  ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("  ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("  ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("=>", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("=>", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("=>", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("=>", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("==", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("==", ::testing::PrintToString(pb));

  pb = {8, 1};
  EXPECT_EQ(" ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ(" ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ(" ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ(" ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ(" ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ(" ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ(" ", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ(" ", ::testing::PrintToString(pb));
  EXPECT_TRUE(++pb);
  EXPECT_EQ("=", ::testing::PrintToString(pb));
  EXPECT_FALSE(++pb);
  EXPECT_EQ("=", ::testing::PrintToString(pb));

  pb = {8, 0};
  EXPECT_EQ(" ", ::testing::PrintToString(pb));
}

}  // namespace test
}  // namespace biosoup
