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
  Logger(const Logger&) = delete;
  Logger& operator=(const Logger&) = delete;

  static Logger& getInstance(std::ostream& out, const LogLevel level) {
    static Logger instance(out, level);
    return instance;
  }

  static Logger& getInstance() {
    return getInstance(std::cout, LogLevel::INFO);
  }

  ~Logger() {
    {
      std::lock_guard<std::mutex> lg(mutex_);
      done_ = true;
    }
    cv_.notify_one();
    if (thread_.joinable()) thread_.join();
  }

  void log(const std::string& msg, const LogLevel level);

 private:
  Logger(std::ostream& out = std::cout, const LogLevel level = LogLevel::INFO)
      : stream_(out), level_(level), done_(false) {
    thread_ = std::thread(&Logger::thread_callback_, this);
  }

  std::ostream& stream_;
  std::deque<std::string> messages_;
  std::condition_variable cv_;
  std::thread thread_;
  std::mutex mutex_;
  LogLevel level_;
  bool done_;

  void thread_callback_();
  void set_stream_(std::ostream& stream);
  std::string get_time_str_();
  std::string level_to_str_(const LogLevel level);
  std::string format_msg_(const std::string& msg, const LogLevel level);
};

#define LOG_FILE(path, log_level)                        \
  do {                                                   \
    static std::ofstream __log_file__(path);             \
    amsim::Logger::getInstance(__log_file__, log_level); \
  } while(0)

#define LOG_STREAM(stream, log_level) \
  amsim::Logger::getInstance(stream, log_level);

#define LOG_INFO(msg) \
  amsim::Logger::getInstance().log(msg, amsim::LogLevel::INFO)
#define LOG_DEBUG(msg) \
  amsim::Logger::getInstance().log(msg, amsim::LogLevel::DEBUG)
#define LOG_WARNING(msg) \
  amsim::Logger::getInstance().log(msg, amsim::LogLevel::WARNING)
#define LOG_ERROR(msg) \
  amsim::Logger::getInstance().log(msg, amsim::LogLevel::ERROR)

}  // namespace amsim

#endif  // AMSIMCPP_LOGGER_H
