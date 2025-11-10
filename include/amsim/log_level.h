#ifndef AMSIMCPP_LOG_LEVEL_H
#define AMSIMCPP_LOG_LEVEL_H

#pragma once

#include <string>

namespace amsim {

enum LogLevel { DEBUG, INFO, WARNING, ERROR, NONE };

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

inline std::ostream& operator<<(std::ostream& os, LogLevel level) {
  return os << to_string(level);
}

}  // namespace amsim

#endif  // AMSIMCPP_LOG_LEVEL_H
