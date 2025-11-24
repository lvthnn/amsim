#ifndef AMSIMCPP_COMPONENTTYPE_H
#define AMSIMCPP_COMPONENTTYPE_H

#include <string>

namespace amsim {

// @brief Enumeration representing phenotype components
//
// The ComponentType enumeration type represents the four different
// phenotype components modelled in the amsimcpp library, namely, (additive)
// genetic, environmental, vertical, and total. The genetic component is
// modelled using the matrix product of a vector of locus effects with the
// standardised genotype matrix of the population, with the environmental
// component being sampled from a multivariate normal distribution. The vertical
// component is a scaled form of the parental phenotypes, and the total
// component is the sum of all these respective components.
enum ComponentType {
  GENETIC = 0,        ///< Genetic component
  ENVIRONMENTAL = 1,  ///< Environmental component
  VERTICAL = 2,       ///< Vertical (nurture) component
  TOTAL = 3           ///< Total (summated) phenotype component
};

/// @brief Postfix-increment a ComponentType object
///
/// Increments the component type, wrapping around to GENETIC if the value
/// passed corresponds to TOTAL.
///
/// @param type The ComponentType instance
/// @return The original value before increment
inline ComponentType operator++(ComponentType& type, int) {
  ComponentType old = type;
  type = (type == ComponentType::TOTAL)
             ? ComponentType::GENETIC
             : ComponentType(static_cast<int>(type) + 1);
  return old;
}

/// @brief Prefix-increment a ComponentType object
///
/// Increments the component type, wrapping around to GENETIC if the value
/// passed corresponds to TOTAL.
///
/// @param type The ComponentType instance
/// @return The incremented ComponentType
inline ComponentType& operator++(ComponentType& type) {
  type = (type == ComponentType::TOTAL)
             ? ComponentType::GENETIC
             : ComponentType(static_cast<int>(type) + 1);
  return type;
}

/// @brief Convert a ComponentType to string
///
/// @param type The ComponentType instance
/// @return The string representation (one of "gen", "env", "vert", or "tot")
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

/// @brief Prints a ComponentType to a stream
///
/// Converts the enum to its string representation and writes
/// it to the given stream.
///
/// @param os The output stream
/// @param type The component type to print
/// @return The same output stream, allowing chaining
inline std::ostream& operator<<(std::ostream& os, ComponentType type) {
  return os << to_string(type);
}

}  // namespace amsim

#endif  // AMSIMCPP_COMPONENTTYPE_H
