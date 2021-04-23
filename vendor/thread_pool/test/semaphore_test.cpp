// Copyright (c) 2020 Robert Vaser

#include "thread_pool/semaphore.hpp"

#include <thread>  // NOLINT

#include "gtest/gtest.h"

namespace thread_pool {
namespace test {

TEST(ThreadPoolSemaphoreTest, Barrier) {
  Semaphore s{0}, b{0};
  auto check = [&] () -> void {
    s.Signal();
    b.Wait();
  };

  std::thread t1{check};
  std::thread t2{check};
  std::thread t3{check};

  s.Wait();
  s.Wait();
  s.Wait();
  // Release barrier
  b.Signal();
  b.Signal();
  b.Signal();

  t1.join();
  t2.join();
  t3.join();
}

}  // namespace test
}  // namespace thread_pool
