#ifndef AMSIMCPP_COMPONENTTYPE_H
#define AMSIMCPP_COMPONENTTYPE_H

namespace amsim {

  enum ComponentType {
    GENETIC       = 0,
    ENVIRONMENTAL = 1,
    VERTICAL      = 2,
    TOTAL         = 3
  };

  inline ComponentType operator++(ComponentType &type, int) {
    type = (type == ComponentType::TOTAL) ? ComponentType::GENETIC : ComponentType(int(type) + 1);
    return type;
  }

}

#endif //AMSIMCPP_COMPONENTTYPE_H
