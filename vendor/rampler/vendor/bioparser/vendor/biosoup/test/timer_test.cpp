// Copyright (c) 2020 Robert Vaser

#include "biosoup/timer.hpp"

#include <thread>  // NOLINT

#include "gtest/gtest.h"

namespace biosoup {
namespace test {

bool ExpectWithinInclusive(double lower, double upper, double value) {
  return lower <= value && value <= upper;
}

TEST(BiosoupTimerTest, Start) {
  Timer t{};
  t.Start();
  std::this_thread::sleep_for(std::chrono::milliseconds(256));
  EXPECT_TRUE(ExpectWithinInclusive(0.256, 0.384, t.Stop()));
  t.Start();
  std::this_thread::sleep_for(std::chrono::milliseconds(256));
  t.Start();
  EXPECT_TRUE(ExpectWithinInclusive(0, 0.001, t.Stop()));
}

TEST(BiosoupTimerTest, Stop) {
  Timer t{};
  EXPECT_TRUE(ExpectWithinInclusive(0, 0.001, t.Stop()));
  EXPECT_EQ(0, t.elapsed_time());
  t.Start();
  std::this_thread::sleep_for(std::chrono::milliseconds(256));
  t.Stop();
  std::this_thread::sleep_for(std::chrono::milliseconds(256));
  EXPECT_TRUE(ExpectWithinInclusive(0.256, 0.384, t.elapsed_time()));
  EXPECT_TRUE(ExpectWithinInclusive(0, 0.001, t.Stop()));
  EXPECT_TRUE(ExpectWithinInclusive(0.256, 0.384, t.elapsed_time()));
  t.Start();
  std::this_thread::sleep_for(std::chrono::milliseconds(256));
  EXPECT_TRUE(ExpectWithinInclusive(0.256, 0.384, t.Stop()));
  EXPECT_TRUE(ExpectWithinInclusive(0.512, 0.768, t.elapsed_time()));
}

TEST(BiosoupTimerTest, Lap) {
  Timer t{};
  t.Start();
  EXPECT_TRUE(ExpectWithinInclusive(0, 0.001, t.Lap()));
  std::this_thread::sleep_for(std::chrono::milliseconds(256));
  EXPECT_TRUE(ExpectWithinInclusive(0.256, 0.384, t.Lap()));
  EXPECT_EQ(0, t.elapsed_time());
  EXPECT_TRUE(ExpectWithinInclusive(0.256, 0.384, t.Stop()));
  EXPECT_TRUE(ExpectWithinInclusive(0, 0.001, t.Lap()));
}

TEST(BiosoupTimerTest, Reset) {
  Timer t{};
  t.Start();
  std::this_thread::sleep_for(std::chrono::milliseconds(256));
  t.Reset();
  t.Stop();
  std::this_thread::sleep_for(std::chrono::milliseconds(256));
  t.Start();
  std::this_thread::sleep_for(std::chrono::milliseconds(256));
  t.Stop();
  EXPECT_TRUE(ExpectWithinInclusive(0.256, 0.384, t.elapsed_time()));
  t.Reset();
  EXPECT_EQ(0, t.elapsed_time());
}

}  // namespace test
}  // namespace biosoup
