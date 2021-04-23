// Copyright (c) 2020 Robert Vaser

#ifndef BIOSOUP_PROGRESS_BAR_HPP_
#define BIOSOUP_PROGRESS_BAR_HPP_

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

namespace biosoup {

class ProgressBar {
 public:
  // num_ticks element of [1, num_events]
  ProgressBar(std::uint32_t num_events, std::uint32_t num_ticks)
      : bar_(std::max(std::min(num_ticks, num_events), 1U), ' '),
        bar_ptr_(0),
        num_events_(num_events),
        event_counter_(0),
        events_per_tick_(num_events / static_cast<double>(bar_.size())) {}

  ProgressBar(const ProgressBar&) = default;
  ProgressBar& operator=(const ProgressBar&) = default;

  ProgressBar(ProgressBar&&) = default;
  ProgressBar& operator=(ProgressBar&&) = default;

  ~ProgressBar() = default;

  std::uint32_t event_counter() const {
    return event_counter_;
  }

  std::uint32_t num_events() const {
    return num_events_;
  }

  bool operator++() {  // true if a tick is added
    ++event_counter_;
    if (event_counter_ > num_events_) {
      return false;
    } else if (event_counter_ == num_events_) {
      while (bar_ptr_ < bar_.size()) {
        bar_[bar_ptr_++] = '=';
      }
      return true;
    } else if (event_counter_ >= (bar_ptr_ + 1) * events_per_tick_) {
      bar_[bar_ptr_++] = '=';
      if (bar_ptr_ < bar_.size()) {
        bar_[bar_ptr_] = '>';
      }
      return true;
    }
    return false;
  }

  // to decrease amount of output, use when ++ returns true
  friend std::ostream& operator<<(std::ostream& os, const ProgressBar& pb) {
    return os << pb.bar_;
  }

 private:
  std::string bar_;
  std::uint32_t bar_ptr_;
  std::uint32_t num_events_;
  std::uint32_t event_counter_;
  double events_per_tick_;
};

}  // namespace biosoup

#endif  // BIOSOUP_PROGRESS_BAR_HPP_
