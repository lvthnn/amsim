#ifndef AMSIMCPP_LOG_LEVEL_H
#define AMSIMCPP_LOG_LEVEL_H

#pragma once

#include <string>

enum LogLevel { DEBUG, INFO, WARNING, ERROR, CRITICAL, NONE };

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
    case LogLevel::CRITICAL:
      return "CRITICAL";
    case LogLevel::NONE:
      return "NONE";
  }
  __builtin_unreachable();
}

inline std::ostream& operator<<(std::ostream& os, LogLevel level) {
  return os << to_string(level);
}

#endif  // AMSIMCPP_LOG_LEVEL_H
