#include <amsim/logger_timer.h>

#include <iomanip>
#include <sstream>

namespace amsim {

LoggerTimer::LoggerTimer(std::string label)
    : label_(std::move(label)), start_(Clock::now()), last_tick_(start_) {}

std::string LoggerTimer::tick(const std::string& message) {
  auto now = Clock::now();
  auto delta = std::chrono::duration<double>(now - last_tick_).count();
  last_tick_ = now;
  std::ostringstream oss;
  if (!label_.empty()) oss << "[" << label_ << "] ";
  oss << message << " (" << std::fixed << std::setprecision(3) << delta
      << " s)";
  return oss.str();
}

}  // namespace amsim
