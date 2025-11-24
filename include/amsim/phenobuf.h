#ifndef AMSIMCPP_PHENOBUF_H
#define AMSIMCPP_PHENOBUF_H

#include <amsim/component_type.h>

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

namespace amsim {
/// @brief Buffer for storing phenotype component values
///
/// PhenoBuf manages storage for phenotype component values (genetic,
/// environmental, vertical, and total) for multiple phenotypes. It optionally
/// supports latent phenotype values for computing correlations. The buffer
/// tracks which phenotype slots are occupied.
class PhenoBuf {
 public:
  /// @brief Construct a new PhenoBuf instance
  ///
  /// @param n_ind Number of individuals in the population
  /// @param n_pheno Number of phenotypes to allocate space for
  /// @param require_lat Whether to allocate latent phenotype buffer
  PhenoBuf(std::size_t n_ind, std::size_t n_pheno, bool require_lat);

  /// @brief Access phenotype component buffer for a specific phenotype
  ///
  /// @param id Phenotype index
  /// @param type Component type to access
  /// @return Const pointer to component value array
  const double* operator()(std::size_t id, ComponentType type) const {
    return &buf_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
  }

  /// @brief Access phenotype component buffer for a specific phenotype
  ///
  /// @param id Phenotype index
  /// @param type Component type to access
  /// @return Pointer to component value array
  double* operator()(std::size_t id, ComponentType type) {
    return &buf_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
  }

  /// @brief Access phenotype component buffer for all phenotypes
  ///
  /// @param type Component type to access
  /// @return Pointer to component value array
  double* operator()(ComponentType type) {
    return &buf_[n_ind_ * n_pheno_ * static_cast<int>(type)];
  }

  /// @brief Access latent phenotype buffer for a specific phenotype
  ///
  /// @param id Phenotype index
  /// @param type Component type to access
  /// @return Const pointer to latent value array
  const double* latent(std::size_t id, ComponentType type) const {
    if (buf_lat_.empty()) throw std::runtime_error("latent buffer not in use");
    return &buf_lat_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
  }

  /// @brief Access latent phenotype buffer for all phenotypes
  ///
  /// @param type Component type to access
  /// @return Pointer to latent value array
  double* latent(ComponentType type) {
    if (buf_lat_.empty()) throw std::runtime_error("latent buffer not in use");
    return &buf_lat_[n_ind_ * n_pheno_ * static_cast<int>(type)];
  }

  /// @brief Return the number of individuals
  /// @return Number of individuals
  std::size_t n_ind() const noexcept { return n_ind_; }

  /// @brief Check if a phenotype slot is occupied
  /// @param id Phenotype index
  /// @return True if slot is occupied
  bool occupied(std::size_t id) const noexcept { return occupied_[id]; }

  /// @brief Check if latent buffer is allocated
  /// @return True if latent buffer exists
  bool has_lat() const noexcept { return !buf_lat_.empty(); }

  /// @brief Find first unoccupied phenotype slot
  /// @return Index of unoccupied slot, or nullopt if all occupied
  std::optional<std::size_t> unoccupied() const;

  /// @brief Mark a phenotype slot as occupied
  /// @param id Phenotype index to mark as occupied
  void occupy(std::size_t);

  /// @brief Compute latent phenotypes from SVD matrices
  ///
  /// @param U Left singular vectors
  /// @param VT Right singular vectors (transposed)
  void score_latent(
      const std::vector<double>& U, const std::vector<double>& VT);

 private:
  const std::size_t n_ind_;      ///< Number of individuals
  const std::size_t n_pheno_;    ///< Number of phenotypes
  std::vector<double> buf_;      ///< Main phenotype component buffer
  std::vector<double> buf_lat_;  ///< Latent phenotype buffer (optional)
  std::vector<bool> occupied_;   ///< Phenotype slot occupancy flags
};

}  // namespace amsim

#endif  // AMSIMCPP_PHENOBUF_H
