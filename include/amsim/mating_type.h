#ifndef AMSIMCPP_MATING_TYPE_H
#define AMSIMCPP_MATING_TYPE_H

namespace amsim {

/// @brief Enumeration of mating model types
///
/// The MatingType enumeration specifies the mating strategy used in the
/// simulation. RANDOM assigns mates randomly, while ASSORTATIVE uses
/// simulated annealing to achieve target phenotypic correlations.
enum class MatingType {
  RANDOM,      ///< Random mating
  ASSORTATIVE  ///< Assortative mating based on phenotypes
};

}  // namespace amsim

#endif  // AMSIMCPP_MATING_TYPE_H
