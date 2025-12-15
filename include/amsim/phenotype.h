#ifndef AMSIMCPP_PHENOTYPE_H
#define AMSIMCPP_PHENOTYPE_H

#include <amsim/component_type.h>
#include <amsim/genome.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/utils.h>

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace amsim {

/// @brief Class representing a phenotypic trait
///
/// Phenotype models a trait as the sum of genetic, environmental, and vertical
/// (nurture) components. It manages trait values for all individuals, computes
/// component statistics, and handles genetic scoring and vertical transmission.
///
/// @see PhenoBuf, PhenoArch
class Phenotype {
 public:
  /// @brief Construct a new Phenotype instance
  ///
  /// @param buf Phenotype buffer for storing component values
  /// @param arch Phenotype architecture defining loci and effects
  /// @param name Name of the phenotype
  /// @param h2_gen Narrow-sense heritability (genetic component)
  /// @param h2_env Environmental variance proportion
  /// @param h2_vert Vertical transmission (nurture) variance proportion
  /// @param rvert_pat Paternal vertical transmission ratio
  /// @param rvert_mat Maternal vertical transmission ratio
  /// @param mate_cor Between-mate correlation on this phenotype
  /// @param id Optional phenotype identifier
  Phenotype(
      PhenoBuf& buf,
      PhenoArch& arch,
      std::string name,
      std::optional<std::size_t> id = std::nullopt);

  /// @brief Return the phenotype name
  /// @return The phenotype name
  std::string name() const noexcept { return name_; }

  /// @brief Return the number of individuals
  /// @return Number of individuals in the population
  std::size_t n_ind() const noexcept { return n_ind_; }

  /// @brief Return the loci contributing to this phenotype
  /// @return Vector of locus indices
  std::vector<std::size_t> loci() const& noexcept { return loci_; }

  /// @brief Return the narrow-sense heritability
  /// @return Genetic variance proportion
  double h2_gen() const noexcept { return h2_gen_; }

  /// @brief Return the vertical transmission proportion
  /// @return Vertical variance proportion
  double h2_vert() const noexcept { return h2_vert_; }

  /// @brief Return the environmental variance proportion
  /// @return Environmental variance proportion
  double h2_env() const noexcept { return h2_env_; }

  double rvert_pat() const noexcept { return rvert_pat_; }

  double rvert_env() const noexcept { return rvert_env_; }

  /// @brief Access phenotype component value for an individual
  ///
  /// @param id Individual index
  /// @param type Component type to access
  /// @return Reference to the component value
  const double& operator()(std::size_t id, ComponentType type) const {
    if (id >= n_ind_)
      throw std::runtime_error("attempting out-of-bounds access of phenotype");
    switch (type) {
      case ComponentType::GENETIC:
        return ptr_gen_[id];
      case ComponentType::ENVIRONMENTAL:
        return ptr_env_[id];
      case ComponentType::VERTICAL:
        return ptr_vert_[id];
      case ComponentType::TOTAL:
        return ptr_tot_[id];
    }
  }

  /// @brief Access total phenotype value for an individual
  ///
  /// @param id Individual index
  /// @return Reference to the total phenotype value
  const double& operator()(std::size_t id) const {
    if (id >= n_ind_)
      throw std::runtime_error("attempting out-of-bounds access of phenotype");
    return ptr_tot_[id];
  }

  /// @brief Access pointer to a phenotype component array
  ///
  /// @param type Component type to access
  /// @return Pointer to the component value array
  const double* operator()(ComponentType type) const {
    switch (type) {
      case ComponentType::GENETIC:
        return ptr_gen_;
      case ComponentType::ENVIRONMENTAL:
        return ptr_env_;
      case ComponentType::VERTICAL:
        return ptr_vert_;
      case ComponentType::TOTAL:
        return ptr_tot_;
    }
    __builtin_unreachable();
  }

  /// @brief Access pointer to a mutable phenotype component array
  ///
  /// @param type Component type to access
  /// @return Pointer to the component value array
  double* operator()(ComponentType type) {
    switch (type) {
      case ComponentType::GENETIC:
        return ptr_gen_;
      case ComponentType::ENVIRONMENTAL:
        return ptr_env_;
      case ComponentType::VERTICAL:
        return ptr_vert_;
      case ComponentType::TOTAL:
        return ptr_tot_;
    }
    __builtin_unreachable();
  }

  /// @brief Return mean of a phenotype component
  ///
  /// @param type Component type
  /// @return Component mean
  double comp_mean(ComponentType type) const { return comp_means_[type]; }

  /// @brief Return variance of a phenotype component
  ///
  /// @param type Component type
  /// @return Component variance
  double comp_var(ComponentType type) const { return comp_vars_[type]; }

  /// @brief Transmit vertical component from parents to offspring
  ///
  /// @param matching Vector of mate pair indices
  void transmit_vert(std::vector<std::size_t> matching);

  /// @brief Compute total phenotype as sum of all components
  void score_tot() {
    for (std::size_t ind = 0; ind < n_ind_; ++ind)
      ptr_tot_[ind] = ptr_gen_[ind] + ptr_env_[ind] + ptr_vert_[ind];
  }

  /// @brief Compute genetic component using bitwise operations
  ///
  /// @param genome Genome instance containing genotype data
  void score_bitwise(Genome& genome);

  /// @brief Compute genetic component using tiled operations
  ///
  /// @param genome Genome instance containing genotype data
  void score_tiled(Genome& genome);

  /// @brief Compute phenotype genetic component
  ///
  /// @param genome Genome instance containing genotype data
  void score(Genome& genome);

  /// @brief Compute component means and variances
  void compute_stats();

 private:
  const std::string name_;               ///< Phenotype name
  const std::size_t id_;                 ///< Phenotype identifier
  const std::size_t n_ind_;              ///< Number of individuals
  const std::vector<std::size_t> loci_;  ///< Loci contributing to trait
  const double* loc_effects_;            ///< Effect sizes for each locus
  const double h2_gen_;                  ///< Narrow-sense heritability
  const double h2_env_;                  ///< Environmental variance proportion
  const double h2_vert_;                 ///< Vertical transmission proportion
  const double rvert_pat_;  ///< Paternal vertical transmission proportion
  const double rvert_env_;  ///< Environmental component transmission proportion

  double* ptr_gen_;   ///< Pointer to genetic component values
  double* ptr_env_;   ///< Pointer to environmental component values
  double* ptr_vert_;  ///< Pointer to vertical component values
  double* ptr_tot_;   ///< Pointer to total phenotype values

  double vert_var_;
  bool vert_lock_ = false;

  std::array<double, 4> comp_means_;  ///< Component means
  std::array<double, 4> comp_vars_;   ///< Component variances
};

/// @brief Type alias for a vector of Phenotype instances
using PhenotypeList = std::vector<Phenotype>;

}  // namespace amsim

#endif  // AMSIMCPP_PHENOTYPE_H
