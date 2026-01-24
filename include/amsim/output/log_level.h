#ifndef AMSIMCPP_LOG_LEVEL_H
#define AMSIMCPP_LOG_LEVEL_H

#include <string>

namespace amsim {

/// @brief Enumeration type representing the level of detail used by logger
enum LogLevel { DEBUG, INFO, WARNING, ERROR, NONE };

/// @brief Convert a LogLevel to string
///
/// @param type The LogLevel instance
/// @return The string representation (one of "DEBUG", "INFO", "WARNING",
///   "ERROR", or "NONE")
inline std::string to_string(LogLevel level) {
  switch (level) {
    case LogLevel::DEBUG:
      return "DEBUG";
    case LogLevel::INFO:
      return "INFO";
    case LogLevel::WARNING:
      return "WARNING";
    case LogLevel::ERROR:
      return "ERROR";
    case LogLevel::NONE:
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

}  // namespace amsim

#endif  // AMSIMCPP_LOG_LEVEL_H
