#ifndef AMSIMCPP_PHENOARCH_H
#define AMSIMCPP_PHENOARCH_H

#include <amsim/rng.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace amsim {

/// @brief Genetic architecture for multiple phenotypes
///
/// PhenoArch manages the genetic architecture of phenotypes, including locus
/// assignments, effect sizes, and genetic/environmental correlations between
/// traits. It optimizes locus assignments to achieve target genetic
/// correlations and generates correlated environmental effects.
class PhenoArch {
 public:
  /// @brief Construct a new PhenoArch instance
  ///
  /// @param n_pheno Number of phenotypes
  /// @param n_loc_total Total number of loci modelled
  /// @param n_loc Number of loci per phenotype
  /// @param h2_gen Narrow-sense heritabilities for each phenotype
  /// @param h2_env Environmental variance proportions
  /// @param gen_cor Target genetic correlation matrix
  /// @param env_cor Target environmental correlation matrix
  /// @param rng A seeded RNG instance
  PhenoArch(
      std::size_t n_pheno,
      std::size_t n_loc_total,
      std::vector<std::size_t> n_loc,
      std::vector<double> h2_gen,
      std::vector<double> h2_env,
      std::vector<double> gen_cor,
      std::vector<double> env_cor,
      const rng::Xoshiro256ss& rng);

  /// @brief Generate correlated environmental effects
  ///
  /// @param ptr_env Pointer to environmental component buffer
  /// @param n_ind Number of individuals
  void gen_env(double* ptr_env, std::size_t n_ind);

  /// @brief Optimize locus assignments to achieve target correlations
  ///
  /// @param max_it Maximum number of optimization iterations
  /// @param eps Convergence tolerance (default: 1e-12)
  void optim_arch(std::size_t max_it, double eps = 1e-12);

  /// @brief Print genetic correlations for intersecting loci
  ///
  /// @param intersect Vector of locus intersection counts
  void print_correlations(const std::vector<std::size_t>& intersect) const;

  /// @brief Return environmental Cholesky decomposition
  /// @return Environmental correlation Cholesky factors
  std::vector<double> env_chol() const noexcept { return env_chol_; }

  /// @brief Return genetic correlations
  /// @return Target genetic correlation matrix
  std::vector<double> gen_cor() const noexcept { return gen_cor_; }

  /// @brief Get locus mask for a specific phenotype
  ///
  /// @param pheno_id Phenotype index
  /// @return Vector of locus indices assigned to this phenotype
  std::vector<std::size_t> pheno_mask(std::size_t pheno_id) const;

 private:
  const std::size_t n_pheno_;             ///< Number of phenotypes
  const std::size_t n_loc_tot_;           ///< Total number of loci modelled
  const std::vector<std::size_t> n_loc_;  ///< Number of loci per phenotype
  const std::vector<double> h2_gen_;      ///< Narrow-sense heritabilities
  const std::vector<double> h2_env_;   ///< Environmental variance proportions
  const std::vector<double> gen_cor_;  ///< Target genetic correlation matrix

  std::vector<double> env_chol_;         ///< Environmental Cholesky factors
  std::vector<std::uint64_t> loc_mask_;  ///< Locus assignment bit masks

  rng::NormalPolar rng_polar_;    ///< RNG for normal variates
  rng::UniformIntRange rng_unf_;  ///< RNG for uniform integers

  /// @brief Initialize random locus assignment masks
  /// @return Vector of locus masks for each phenotype
  std::vector<std::uint64_t> initMask();

  /// @brief Compute initial energy for optimization
  /// @return Vector of energy values
  std::vector<double> initEnergy() const;

  /// @brief Compute locus effect weights
  /// @return Vector of effect size weights
  std::vector<double> initWeights() const;

  /// @brief Count locus intersections between phenotypes
  ///
  /// @param mask Vector of locus masks
  /// @return Vector of pairwise locus intersection counts
  std::vector<std::size_t> initIntersect(
      const std::vector<std::uint64_t>& mask) const;
};

}  // namespace amsim

#endif  // AMSIMCPP_PHENOARCH_H
