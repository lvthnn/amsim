#pragma once

#include <condition_variable>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>

#include <boost/algorithm/string.hpp>

namespace amsim {

enum LogLevel { Debug, Info, Warning, Error, None };

inline std::string LogLevel_to_string(LogLevel level) {
  switch (level) {
    case LogLevel::Debug:
      return "DEBUG";
    case LogLevel::Info:
      return "INFO";
    case LogLevel::Warning:
      return "WARNING";
    case LogLevel::Error:
      return "ERROR";
    case LogLevel::None:
      return "NONE";
  }
  __builtin_unreachable();
}

inline LogLevel LogLevel_from_string(const std::string& level) {
  std::string upper = boost::to_upper_copy(level);
  if (upper == "DEBUG") return LogLevel::Debug;
  if (upper == "INFO") return LogLevel::Info;
  if (upper == "WARNING") return LogLevel::Warning;
  if (upper == "ERROR") return LogLevel::Error;
  if (upper == "NONE") return LogLevel::None;
  throw std::runtime_error("Unknown LogLevel " + upper);
}

inline std::ostream& operator<<(std::ostream& os, LogLevel level) {
  return os << LogLevel_to_string(level);
}

class Log {
 public:
  Log(const Log&) = delete;
  Log& operator=(const Log&) = delete;

  static Log& get_instance(std::ostream& out, const LogLevel level) {
    static Log instance(out, level);
    return instance;
  }

  static Log& get_instance() { return get_instance(std::cout, LogLevel::Info); }

  ~Log() {
    {
      std::lock_guard<std::mutex> lg(mutex_);
      done_ = true;
    }
    cv_.notify_one();
    if (thread_.joinable()) thread_.join();
  }

  void log(const std::string& msg, LogLevel level);
  static void file(const std::filesystem::path& path, LogLevel level);
  static void stream(std::ostream& stream, LogLevel level);
  static void debug(const std::string& msg);
  static void info(const std::string& msg);
  static void warning(const std::string& msg);
  static void error(const std::string& msg);

 private:
  explicit Log(
      std::ostream& out = std::cout, const LogLevel level = LogLevel::Info)
      : stream_(out), level_(level) {
    thread_ = std::thread(&Log::threadCallback, this);
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

inline void Log::threadCallback() {
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

inline std::string Log::getTimeStr() {
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

inline std::string Log::formatMsg(
    const std::string& msg, const LogLevel level) {
  std::string time_str = getTimeStr();
  std::ostringstream msg_format;

  msg_format << "[" << LogLevel_to_string(level) << "] "
             << "[" << time_str << "] " << "[thread "
             << std::this_thread::get_id() << "] " << msg;

  return msg_format.str();
}

inline void Log::log(const std::string& msg, const LogLevel level) {
  if (level < level_) return;
  {
    std::lock_guard<std::mutex> lg(mutex_);
    std::string msg_format = formatMsg(msg, level);
    messages_.push_back(msg_format);
  }
  cv_.notify_one();
}

inline void Log::file(const std::filesystem::path& path, LogLevel log_level) {
  static std::fstream log_file(path, std::ios::out);
  Log::get_instance(log_file, log_level);
}

inline void Log::stream(std::ostream& stream, LogLevel log_level) {
  Log::get_instance(stream, log_level);
}

inline void Log::debug(const std::string& msg) {
  Log::get_instance().log(msg, LogLevel::Debug);
}

inline void Log::info(const std::string& msg) {
  Log::get_instance().log(msg, LogLevel::Info);
}

inline void Log::warning(const std::string& msg) {
  Log::get_instance().log(msg, LogLevel::Warning);
}

inline void Log::error(const std::string& msg) {
  Log::get_instance().log(msg, LogLevel::Error);
}

}  // namespace amsim
