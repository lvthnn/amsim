#ifndef AMSIMCPP_COMPONENTTYPE_H
#define AMSIMCPP_COMPONENTTYPE_H

#pragma once

#include <string>

namespace amsim {

enum ComponentType { GENETIC = 0, ENVIRONMENTAL = 1, VERTICAL = 2, TOTAL = 3 };

inline ComponentType operator++(ComponentType& type, int) {
  ComponentType old = type;
  type = (type == ComponentType::TOTAL) ? ComponentType::GENETIC
                                        : ComponentType(int(type) + 1);
  return old;
}

inline ComponentType& operator++(ComponentType& type) {
  type = (type == ComponentType::TOTAL) ? ComponentType::GENETIC
                                        : ComponentType(int(type) + 1);
  return type;
}

inline std::string to_string(ComponentType type) {
  switch (type) {
    case ComponentType::GENETIC:
      return "gen";
    case ComponentType::ENVIRONMENTAL:
      return "env";
    case ComponentType::VERTICAL:
      return "vert";
    case ComponentType::TOTAL:
      return "tot";
  }
  __builtin_unreachable();
}

inline std::ostream& operator<<(std::ostream& os, ComponentType type) {
  return os << to_string(type);
}

}  // namespace amsim

#endif  // AMSIMCPP_COMPONENTTYPE_H
