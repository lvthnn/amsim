#pragma once

#include <condition_variable>
#include <deque>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>

namespace amsim {

enum LogLevel { Debug, Info, Warning, Error, None };

inline std::string to_string(LogLevel level) {
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

/// @brief Prints a LogLevel to a stream
///
/// Converts the enum to its string representation and writes
/// it to the given stream.
///
/// @param os The output stream
/// @param type The log level to print
/// @return The same output stream, allowing chaining
inline std::ostream& operator<<(std::ostream& os, LogLevel level) {
  return os << to_string(level);
}

/// @brief Timer for logging elapsed time intervals
///
/// LoggerTimer measures and logs time intervals between checkpoints, useful
/// for performance profiling and progress monitoring.
class LoggerTimer {
 public:
  using Clock = std::chrono::high_resolution_clock;  ///< Shorten
  using TimePoint = Clock::time_point;               ///< Time point type

  /// @brief Construct a LoggerTimer
  /// @param label Optional label for timer
  explicit LoggerTimer(std::string label = "")
      : label_(std::move(label)), start_(Clock::now()), last_tick_(start_) {}

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

inline std::string LoggerTimer::tick(const std::string& message) {
  auto now = Clock::now();
  auto delta = std::chrono::duration<double>(now - last_tick_).count();
  last_tick_ = now;
  std::ostringstream oss;
  if (!label_.empty()) oss << "[" << label_ << "] ";
  oss << message << " (" << std::fixed << std::setprecision(3) << delta
      << " s)";
  return oss.str();
}

/// @brief Thread-safe singleton logger
///
/// Logger provides asynchronous, thread-safe logging with configurable
/// verbosity levels. It uses a background thread to write log messages.
class Logger {
 public:
  Logger(const Logger&) = delete;
  Logger& operator=(const Logger&) = delete;

  /// @brief Get Logger instance with custom stream and level
  /// @param out Output stream
  /// @param level Logging level
  /// @return Reference to singleton Logger instance
  static Logger& get_instance(std::ostream& out, const LogLevel level) {
    static Logger instance(out, level);
    return instance;
  }

  /// @brief Get Logger instance with default settings
  /// @return Reference to singleton Logger instance
  static Logger& get_instance() {
    return get_instance(std::cout, LogLevel::Info);
  }

  /// @brief Destructor - flushes messages and joins logger thread
  ~Logger() {
    {
      std::lock_guard<std::mutex> lg(mutex_);
      done_ = true;
    }
    cv_.notify_one();
    if (thread_.joinable()) thread_.join();
  }

  /// @brief Log a message at specified level
  /// @param msg Message to log
  /// @param level Log level for this message
  void log(const std::string& msg, LogLevel level);

 private:
  /// @brief Private constructor for singleton
  /// @param out Output stream
  /// @param level Minimum log level
  explicit Logger(
      std::ostream& out = std::cout, const LogLevel level = LogLevel::Info)
      : stream_(out), level_(level) {
    thread_ = std::thread(&Logger::threadCallback, this);
  }

  std::ostream& stream_;              ///< Output stream
  std::deque<std::string> messages_;  ///< Message queue
  std::condition_variable cv_;  ///< Condition variable for synchronization
  std::thread thread_;          ///< Background logging thread
  std::mutex mutex_;            ///< Mutex for thread safety
  LogLevel level_;              ///< Minimum log level
  bool done_ = false;           ///< Shutdown flag

  /// @brief Background thread callback for writing log messages
  void threadCallback();

  /// @brief Set output stream
  /// @param stream Output stream
  void setStream(std::ostream& stream);

  /// @brief Convert log level to string
  /// @param level Log level
  /// @return String representation
  std::string levelToStr(LogLevel level);

  /// @brief Get current timestamp string
  /// @return Formatted timestamp
  static std::string getTimeStr();

  /// @brief Format log message with timestamp and level
  /// @param msg Message text
  /// @param level Log level
  /// @return Formatted message
  static std::string formatMsg(const std::string& msg, LogLevel level);
};

inline void Logger::threadCallback() {
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

inline std::string Logger::getTimeStr() {
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

inline std::string Logger::formatMsg(
    const std::string& msg, const LogLevel level) {
  std::string time_str = getTimeStr();
  std::ostringstream msg_format;

  msg_format << "[" << to_string(level) << "] "
             << "[" << time_str << "] " << "[thread "
             << std::this_thread::get_id() << "] " << msg;

  return msg_format.str();
}

inline void Logger::log(const std::string& msg, const LogLevel level) {
  if (level < level_) return;
  {
    std::lock_guard<std::mutex> lg(mutex_);
    std::string msg_format = formatMsg(msg, level);
    messages_.push_back(msg_format);
  }
  cv_.notify_one();
}

/// @brief Initialize logger to write to a file
/// @param path File path
/// @param log_level Minimum log level
#define LOG_FILE(path, log_level)                         \
  do {                                                    \
    static std::ofstream __log_file__(path);              \
    amsim::Logger::get_instance(__log_file__, log_level); \
  } while (false)

/// @brief Initialize logger to write to a stream
/// @param stream Output stream
/// @param log_level Minimum log level
#define LOG_STREAM(stream, log_level) \
  amsim::Logger::get_instance(stream, log_level);

/// @brief Log an informational message
/// @param msg Message to log
#define LOG_INFO(msg) \
  amsim::Logger::get_instance().log(msg, amsim::LogLevel::Info)

/// @brief Log a debug message
/// @param msg Message to log
#define LOG_DEBUG(msg) \
  amsim::Logger::get_instance().log(msg, amsim::LogLevel::Debug)

/// @brief Log a warning message
/// @param msg Message to log
#define LOG_WARNING(msg) \
  amsim::Logger::get_instance().log(msg, amsim::LogLevel::Warning)

/// @brief Log an error message
/// @param msg Message to log
#define LOG_ERROR(msg) \
  amsim::Logger::get_instance().log(msg, amsim::LogLevel::Error)

}  // namespace amsim
