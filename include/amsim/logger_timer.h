#ifndef AMSIMCPP_LOGGER_TIMER_H
#define AMSIMCPP_LOGGER_TIMER_H

#include <chrono>
#include <string>

namespace amsim {

class LoggerTimer {
 public:
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = Clock::time_point;

  explicit LoggerTimer(std::string label = "");
  std::string tick(const std::string& message = "");
  ~LoggerTimer() = default;

 private:
  std::string label_;
  TimePoint start_;
  TimePoint last_tick_;
};

}  // namespace amsim

#endif  // AMSIMCPP_LOGGER_TIMER_H
