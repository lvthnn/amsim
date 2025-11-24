#ifndef AMSIMCPP_LOGGER_TIMER_H
#define AMSIMCPP_LOGGER_TIMER_H

#include <chrono>
#include <string>

namespace amsim {

/// @brief Timer for logging elapsed time intervals
///
/// LoggerTimer measures and logs time intervals between checkpoints, useful
/// for performance profiling and progress monitoring.
class LoggerTimer {
 public:
  using Clock = std::chrono::high_resolution_clock;  ///< Clock type
  using TimePoint = Clock::time_point;               ///< Time point type

  /// @brief Construct a LoggerTimer
  /// @param label Optional label for timer
  explicit LoggerTimer(std::string label = "");

  /// @brief Record a checkpoint and return elapsed time message
  /// @param message Optional message for this checkpoint
  /// @return Formatted string with elapsed time since last tick
  std::string tick(const std::string& message = "");

  ~LoggerTimer() = default;

 private:
  std::string label_;    ///< Timer label
  TimePoint start_;      ///< Start time
  TimePoint last_tick_;  ///< Last checkpoint time
};

}  // namespace amsim

#endif  // AMSIMCPP_LOGGER_TIMER_H
