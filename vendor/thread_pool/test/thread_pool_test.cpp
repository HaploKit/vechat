// Copyright (c) 2020 Robert Vaser

#include "thread_pool/thread_pool.hpp"

#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "gtest/gtest.h"

namespace thread_pool {
namespace test {

TEST(ThreadPoolThreadPoolTest, Submit) {
  std::function<std::uint32_t(std::uint32_t)> fibonacci =
      [&fibonacci] (std::uint32_t i) -> std::uint32_t {
    if (i == 1 || i == 2) {
      return 1;
    }
    return fibonacci(i - 1) + fibonacci(i - 2);
  };

  ThreadPool tp{};

  std::vector<std::future<std::uint32_t>> f;
  for (std::uint32_t i = 0; i < tp.num_threads(); ++i) {
    f.emplace_back(tp.Submit(fibonacci, 42));
  }
  for (auto& it : f) {
    EXPECT_EQ(267914296, it.get());
  }
}

TEST(ThreadPoolThreadPoolTest, ThreadIds) {
  ThreadPool tp{};
  EXPECT_EQ(tp.num_threads(), tp.thread_ids().size());

  Semaphore s{0}, b{0};
  auto check = [&] () -> std::uint32_t {
    EXPECT_EQ(1, tp.thread_ids().count(std::this_thread::get_id()));
    s.Signal();
    b.Wait();
    return tp.thread_ids().at(std::this_thread::get_id());
  };

  std::vector<std::future<std::uint32_t>> f;
  for (std::uint32_t i = 0; i < tp.num_threads(); ++i) {
    f.emplace_back(tp.Submit(check));
  }

  for (std::uint32_t i = 0; i < tp.num_threads(); ++i) {
    s.Wait();
  }
  for (std::uint32_t i = 0; i < tp.num_threads(); ++i) {
    b.Signal();
  }

  std::unordered_set<std::uint32_t> ts;
  for (auto& it : f) {
    ts.emplace(it.get());
  }
  EXPECT_EQ(tp.num_threads(), ts.size());
}

}  // namespace test
}  // namespace thread_pool
