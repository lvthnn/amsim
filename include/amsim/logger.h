#ifndef AMSIMCPP_LOGGER_H
#define AMSIMCPP_LOGGER_H

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

  void log(const std::string& msg, LogLevel level);

 private:
  explicit Logger(
      std::ostream& out = std::cout, const LogLevel level = LogLevel::INFO)
      : stream_(out), level_(level) {
    thread_ = std::thread(&Logger::threadCallback, this);
  }

  std::ostream& stream_;
  std::deque<std::string> messages_;
  std::condition_variable cv_;
  std::thread thread_;
  std::mutex mutex_;
  LogLevel level_;
  bool done_ = false;

  void threadCallback();
  void setStream(std::ostream& stream);
  std::string levelToStr(LogLevel level);
  static std::string getTimeStr();
  static std::string formatMsg(const std::string& msg, LogLevel level);
};

#define LOG_FILE(path, log_level)                        \
  do {                                                   \
    static std::ofstream __log_file__(path);             \
    amsim::Logger::getInstance(__log_file__, log_level); \
  } while (false)

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
