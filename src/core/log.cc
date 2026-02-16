#include <amsim/core/log.h>

#include <chrono>
#include <condition_variable>
#include <deque>
#include <iomanip>
#include <iostream>
#include <mutex>
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

void Logger::threadCallback() {
  std::unique_lock<std::mutex> lk(mutex_);
  while (!done_ || !messages_.empty()) {
    cv_.wait(lk, [this] { return done_ || !messages_.empty(); });

    while (!messages_.empty()) {
      std::string msg = std::move(messages_.front());
      messages_.pop_front();

      lk.unlock();
      stream_ << msg << "\n";
      lk.lock();
    }
  }
}

std::string Logger::getTimeStr() {
  const auto now = std::chrono::system_clock::now();
  const auto tse = now.time_since_epoch();
  const auto ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(tse) % 1000;
  std::time_t t = std::chrono::system_clock::to_time_t(now);
  std::tm tm_now;

#if defined(_WIN32) || defined(_WIN64)
  localtime_s(&tm_now, &t);
#else
  localtime_r(&t, &tm_now);
#endif

  std::ostringstream oss;

  oss << std::put_time(&tm_now, "%Y-%m-%d %H:%M:%S") << '.' << std::setfill('0')
      << std::setw(3) << ms.count();

  return oss.str();
}

std::string Logger::formatMsg(const std::string& msg, const LogLevel level) {
  std::string time_str = getTimeStr();
  std::ostringstream msg_format;

  msg_format << "[" << to_string(level) << "] "
             << "[" << time_str << "] " << "[thread "
             << std::this_thread::get_id() << "] " << msg;

  return msg_format.str();
}

void Logger::log(const std::string& msg, const LogLevel level) {
  if (level < level_) return;
  {
    std::lock_guard<std::mutex> lg(mutex_);
    std::string msg_format = formatMsg(msg, level);
    messages_.push_back(msg_format);
  }
  cv_.notify_one();
}

}  // namespace amsim

