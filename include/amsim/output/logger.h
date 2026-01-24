#ifndef AMSIMCPP_LOGGER_H
#define AMSIMCPP_LOGGER_H

#include <amsim/log_level.h>

#include <condition_variable>
#include <deque>
#include <iostream>
#include <mutex>
#include <thread>

namespace amsim {

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
    return get_instance(std::cout, LogLevel::INFO);
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
      std::ostream& out = std::cout, const LogLevel level = LogLevel::INFO)
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
  amsim::Logger::get_instance().log(msg, amsim::LogLevel::INFO)

/// @brief Log a debug message
/// @param msg Message to log
#define LOG_DEBUG(msg) \
  amsim::Logger::get_instance().log(msg, amsim::LogLevel::DEBUG)

/// @brief Log a warning message
/// @param msg Message to log
#define LOG_WARNING(msg) \
  amsim::Logger::get_instance().log(msg, amsim::LogLevel::WARNING)

/// @brief Log an error message
/// @param msg Message to log
#define LOG_ERROR(msg) \
  amsim::Logger::get_instance().log(msg, amsim::LogLevel::ERROR)

}  // namespace amsim

#endif  // AMSIMCPP_LOGGER_H
