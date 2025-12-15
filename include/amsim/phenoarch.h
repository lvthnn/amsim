#ifndef AMSIMCPP_PHENOARCH_H
#define AMSIMCPP_PHENOARCH_H

#include <amsim/logger.h>
#include <amsim/rng.h>
#include <amsim/utils.h>

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
  /// @param h2_vert Vertical variance proportions
  /// @param gen_cor Target genetic correlation matrix
  /// @param env_cor Target environmental correlation matrix
  PhenoArch(
      std::size_t n_pheno,
      std::size_t n_loc_total,
      std::vector<std::size_t> n_loc,
      std::vector<std::string> pheno_names,
      std::vector<double> v_h2_gen,
      std::vector<double> v_h2_env,
      std::vector<double> v_h2_vert,
      std::vector<double> v_rvert_pat,
      std::vector<double> v_rvert_env,
      std::vector<double> gen_cor,
      std::vector<double> env_cor);

  /// @brief Generate correlated environmental effects
  ///
  /// @param ptr_env Pointer to environmental component buffer
  /// @param n_ind Number of individuals in population
  void gen_env(double* ptr_env, std::size_t n_ind);

  /// @brief Generate initial vertical components
  ///
  /// @param ptr_vert Pointer to phenotype vertical buffer
  /// @param n_ind Number of individuals in population
  /// @param h2_vert Vertical component variance proportion
  static void gen_vert(double* ptr_vert, std::size_t n_ind, double h2_vert);

  /// @brief Optimize locus assignments to achieve target correlations
  ///
  /// @param max_it Maximum number of optimization iterations
  /// @param eps Convergence tolerance (default: 1e-12)
  void optim_arch(std::size_t max_it, double eps = 1e-12);

  /// @brief Return number of phenotypes
  /// @return Number of phenotypes
  std::size_t n_pheno() const { return n_pheno_; }

  /// @brief Return number of causal loci for a phenotype
  ///
  /// @param pheno_id Integer ID of phenotype
  /// @return Number of causal loci for indicated phenotype
  std::size_t n_loc(std::size_t pheno_id) const { return v_n_loc_[pheno_id]; }

  /// @brief Return genetic component variance (heritability) for a phenotype
  ///
  /// @param pheno_id Integer ID of phenotype
  /// @return Genetic component variance for indicated phenotype
  double h2_gen(std::size_t pheno_id) const { return v_h2_gen_[pheno_id]; }

  /// @brief Return environmental component variance for a phenotype
  ///
  /// @param pheno_id Integer ID of phenotype
  /// @return Environmental component variance for indicated phenotype
  double h2_env(std::size_t pheno_id) const { return v_h2_env_[pheno_id]; }

  /// @brief Return vertical component variance for a phenotype
  ///
  /// @param pheno_id Integer ID of phenotype
  /// @return Vertical component variance for indicated phenotype
  double h2_vert(std::size_t pheno_id) const { return v_h2_vert_[pheno_id]; }

  /// @brief Return ratio of paternal vertical transmission
  ///
  /// @param pheno_id Integer ID of phenotype
  /// @return Ratio of paternal vertical transmission for indicated phenotype
  double rvert_pat(std::size_t pheno_id) const {
    return v_rvert_pat_[pheno_id];
  }

  void calibrate_scale(std::size_t pheno_id, double var_vert) {
    if (vert_lock_[pheno_id])
      throw std::runtime_error("noise variance is already calibrated");

    vert_scale_[pheno_id] = std::sqrt(v_h2_vert_[pheno_id] / var_vert);
    vert_lock_[pheno_id] = true;

    LOG_DEBUG(
        "scaling factor for pheno " + std::to_string(pheno_id) + ": " +
        std::to_string(vert_scale_[pheno_id]));
  }

  bool vert_lock(std::size_t pheno_id) const noexcept {
    return vert_lock_[pheno_id];
  }

  double vert_scale(std::size_t pheno_id) const noexcept {
    return vert_scale_[pheno_id];
  }

  /// @brief Return non-transmitted / environmental vertical transmission ratio
  ///
  /// @param pheno_id Integer ID of phenotype
  /// @return Non-transmitted / environmental ratio of vertical transmission
  double rvert_env(std::size_t pheno_id) const {
    return v_rvert_env_[pheno_id];
  }

  /// @brief Return environmental Cholesky decomposition
  /// @return Environmental correlation Cholesky factors
  std::vector<double> env_chol() const noexcept { return env_chol_; }

  /// @brief Return genetic correlations
  /// @return Target genetic correlation matrix
  std::vector<double> gen_cor() const noexcept { return gen_cor_; }

  /// @brief Get locus mask for a phenotype
  ///
  /// @param pheno_id Phenotype index
  /// @return Vector of locus indices assigned to the indicated phenotype
  std::vector<std::size_t> pheno_loc(std::size_t pheno_id) const;

  /// @brief Get the locus effect vector for a phenotype
  ///
  /// @param pheno_id Phenotype index
  /// @return Vector of locus effect sizes for the indicated phenotype
  const double* pheno_effects(std::size_t pheno_id) const;

 private:
  const std::size_t n_pheno_;               ///< Number of phenotypes
  const std::size_t n_loc_tot_;             ///< Total number of loci modelled
  const std::vector<std::size_t> v_n_loc_;  ///< Number of loci per phenotype
  const std::vector<std::string> v_names_;  ///< Phenotype names
  const std::vector<double> v_h2_gen_;      ///< Narrow-sense heritabilities
  const std::vector<double> v_h2_env_;     ///< Environmental component variance
  const std::vector<double> v_h2_vert_;    ///< Vertical component variance
  const std::vector<double> v_rvert_pat_;  ///< Paternal vertical proportion
  const std::vector<double> v_rvert_env_;  ///< Environment vertical proportion
  const std::vector<double> gen_cor_;  ///< Target genetic correlation matrix

  std::vector<double> env_chol_;         ///< Environmental Cholesky factors
  std::vector<std::uint64_t> loc_mask_;  ///< Locus assignment bit masks
  std::vector<double> loc_effects_;      ///< Buffer of phenotype effect sizes

  std::vector<bool> vert_lock_;
  std::vector<double> vert_scale_;

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
