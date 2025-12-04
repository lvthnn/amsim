#ifndef AMSIMCPP_MATING_H
#define AMSIMCPP_MATING_H

#include <amsim/mating_type.h>
#include <amsim/phenotype.h>
#include <amsim/rng.h>

#include <cstddef>
#include <random>
#include <vector>

namespace amsim {
/// @brief Abstract base class for mating models
///
/// MatingModel provides the interface for mating strategies used to determine
/// how individuals pair for reproduction in population genetics simulations.
/// Derived classes implement specific mating strategies by overriding the
/// match() method to produce a vector of mate pairings.
class MatingModel {
 public:
  /// @brief Construct a new MatingModel instance
  ///
  /// @param type The mating type strategy used
  /// @param n_sex Number of individuals in each sex class
  explicit MatingModel(MatingType type, std::size_t n_sex)
      : type_(type), g_(std::random_device{}()), n_sex_(n_sex) {}

  virtual ~MatingModel() = default;

  /// @brief Generate mate pairings
  ///
  /// @return Vector of mate pair indices
  virtual std::vector<std::size_t> match() = 0;

 protected:
  MatingType type_;    ///< Mating type strategy
  std::mt19937 g_;     ///< Random number generator
  std::size_t n_sex_;  ///< Number of individuals in each sex class

  /// @brief Generate a random permutation of mate indices
  ///
  /// @return Vector containing a random permutation of indices from 0 to
  /// n_sex_-1
  std::vector<std::size_t> randState();
};

/// @brief Assortative mating model using simulated annealing
///
/// AssortativeModel implements assortative mating by finding mate pairings
/// that achieve target phenotypic correlations between partners. It uses
/// simulated annealing to optimize the matching arrangement, minimizing the
/// difference between observed and target cross-sex phenotypic correlations.
///
/// @see MatingModel
class AssortativeModel : public MatingModel {
 public:
  /// @brief Construct a new AssortativeModel instance
  ///
  /// @param phenotypes List of phenotypes to use for mate matching
  /// @param cor Vector of target cross-sex phenotypic correlations
  /// @param n_sex Number of individuals in each sex class
  /// @param n_itr Number of iterations for simulated annealing (default: 2e6)
  /// @param temp_init Initial temperature for annealing (default: 0.50)
  /// @param temp_decay Temperature decay factor per iteration (default: 0.9999)
  /// @param tol_inf Tolerance threshold for infinite values (default: 1e-7)
  AssortativeModel(
      const PhenotypeList& phenotypes,
      std::vector<double> cor,
      std::size_t n_sex,
      std::size_t n_itr = 2e6,
      double temp_init = 0.50,
      double temp_decay = 0.9999,
      double tol_inf = 1e-7);

  /// @brief Initialize mate pairing state to a random permutation
  void init_state();

  /// @brief Generate mate pairings using simulated annealing
  ///
  /// @return Vector of mate pair indices optimized for target correlations
  std::vector<std::size_t> match() override;

  /// @brief Update phenotype pointers and recompute statistics
  ///
  /// @param phenotypes Updated list of phenotypes
  void update(const PhenotypeList& phenotypes);

 private:
  /// Pointers to total phenotype values for each phenotype
  std::vector<const double*> ptr_tot_;

  const std::vector<double> cor_;  ///< Target cross-sex phenotypic correlations
  const std::size_t n_pheno_;      ///< Number of phenotypes being matched
  const std::size_t n_sex_;        ///< Number of individuals in each sex class
  const double tol_inf_;           ///< Tolerance threshold for infinite values
  const std::size_t n_itr_;        ///< Number of annealing iterations
  const double temp_init_;         ///< Initial annealing temperature
  const double temp_decay_;        ///< Temperature decay factor

  std::vector<double> male_;    ///< Male phenotype values (flattened)
  std::vector<double> female_;  ///< Female phenotype values (flattened)

  /// @brief Arrange phenotype values based on current mating state
  void arrange();

  /// @brief Compute current cross-sex phenotypic correlations
  ///
  /// @return Vector of correlations in the current mating arrangement
  std::vector<double> computeCor();

  /// @brief Compute change in correlations from swapping two pairings
  ///
  /// @param i0 Index of first individual to swap
  /// @param i1 Index of second individual to swap
  /// @return Vector of correlation changes resulting from the swap
  std::vector<double> computeDelta(std::size_t i0, std::size_t i1);

  /// @brief Compute energy difference for annealing acceptance criterion
  ///
  /// @param cur Current correlation values
  /// @param target Target correlation values
  /// @param delta Change in correlation values
  /// @return Energy difference value
  double computeDiffEnergy(
      const std::vector<double>& cur,
      const std::vector<double>& target,
      const std::vector<double>& delta) const;

 public:
  std::vector<double> cor_S;       ///< Singular values from SVD
  std::vector<double> cor_U;       ///< Left singular vectors from SVD
  std::vector<double> cor_VT;      ///< Right singular vectors from SVD
  std::vector<std::size_t> state;  ///< Current mate pairing state
};
}  // namespace amsim
#endif  // AMSIMCPP_MATING_H
