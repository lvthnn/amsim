#ifndef AMSIMCPP_LOGGER_H
#define AMSIMCPP_LOGGER_H

#pragma once

#include <amsim/log_level.h>

#include <iostream>
#include <mutex>
#include <thread>

namespace amsim {

class Logger {
 public:
  explicit Logger(
      LogLevel level = LogLevel::INFO, std::ostream& out = std::cout)
      : level_(level), out_(out), done_(false) {
    thread_ = std::thread(&Logger::thread_callback_, this);
  }

  ~Logger() {
    {
      std::lock_guard<std::mutex> lg(mutex_);
      done_ = true;
    }
    cv_.notify_one();
    if (thread_.joinable()) thread_.join();
  }

  inline void debug(const std::string& msg) {
    return log_(msg, LogLevel::DEBUG);
  }
  inline void info(const std::string& msg) { return log_(msg, LogLevel::INFO); }
  inline void warning(const std::string& msg) {
    return log_(msg, LogLevel::WARNING);
  }
  inline void error(const std::string& msg) {
    return log_(msg, LogLevel::ERROR);
  }
  inline void critical(const std::string& msg) {
    return log_(msg, LogLevel::CRITICAL);
  }

 private:
  LogLevel level_;
  std::ostream& out_;
  std::deque<std::string> messages_;
  std::condition_variable cv_;
  std::thread thread_;
  std::mutex mutex_;
  bool done_;

  std::string get_time_str_();
  std::string level_to_str_(const LogLevel level);
  std::string format_msg_(const std::string& msg, const LogLevel level);
  void log_(const std::string& msg, const LogLevel level);

  void thread_callback_();
};

}  // namespace amsim

#endif  // AMSIMCPP_LOGGER_H
