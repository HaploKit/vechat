// Copyright (c) 2020 Robert Vaser

#ifndef BIOSOUP_TIMER_HPP_
#define BIOSOUP_TIMER_HPP_

#include <chrono>  // NOLINT
#include <cstdint>

namespace biosoup {

class Timer {
 public:
  Timer()
      : checkpoint_(), elapsed_time_(0) {}

  Timer(const Timer&) = default;
  Timer& operator=(const Timer&) = default;

  Timer(Timer&&) = default;
  Timer& operator=(Timer&&) = default;

  ~Timer() = default;

  double elapsed_time(void) const {
    return elapsed_time_;
  }

  void Start() {
    checkpoint_ = std::chrono::steady_clock::now();
  }

  double Stop() {
    if (checkpoint_.time_since_epoch().count()) {  // Start() was called
      auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - checkpoint_).count();
      checkpoint_ = {};
      elapsed_time_ += duration;
      return duration;
    }
    return 0;
  }

  double Lap() const {
    if (checkpoint_.time_since_epoch().count()) {  // Start() was called
      return std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - checkpoint_).count();
    }
    return 0;
  }

  void Reset() {
    checkpoint_ = {};
    elapsed_time_ = 0;
  }

 private:
  std::chrono::time_point<std::chrono::steady_clock> checkpoint_;
  double elapsed_time_;
};

}  // namespace biosoup

#endif  // BIOSOUP_TIMER_HPP_
