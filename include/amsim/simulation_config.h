#ifndef AMSIMCPP_SIMULATION_CONFIG_H
#define AMSIMCPP_SIMULATION_CONFIG_H

#include <amsim/genome.h>
#include <amsim/haplobuf.h>
#include <amsim/mating.h>
#include <amsim/metric.h>
#include <amsim/metricspec.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/phenotype.h>
#include <amsim/rng.h>

#include <filesystem>
#include <optional>

namespace amsim {

/// @brief Configuration struct for assortative mating simulations
///
/// SimulationConfig uses a builder pattern to configure simulation parameters
/// including genome architecture, phenotypes, mating models, and metrics.
struct SimulationConfig {
  SimulationConfig() = default;

  /// @brief Set simulation parameters
  ///
  /// @param n_gen_ Number of generations to simulate
  /// @param n_ind_ Population size
  /// @param out_dir_ Output directory path
  /// @param rng_seed_ Optional RNG seed
  /// @return Reference to this for method chaining
  SimulationConfig& simulation(
      std::size_t n_gen_,
      std::size_t n_ind_,
      const std::string& out_dir_,
      std::optional<std::uint64_t> rng_seed_);

  /// @brief Set genome parameters
  ///
  /// @param n_loc_ Number of loci
  /// @param v_maf_ Minor allele frequency per locus
  /// @param v_rec_ Recombination probability per locus
  /// @param v_mut_ Mutation probability per locus
  /// @return Reference to this for method chaining
  SimulationConfig& genome(
      std::size_t n_loc_,
      std::vector<double> v_maf_,
      std::vector<double> v_rec_,
      std::vector<double> v_mut_);

  /// @brief Set phenotype parameters
  ///
  /// @param n_pheno_ Number of phenotypes
  /// @param v_name_ Phenotype names
  /// @param v_n_loc_ Loci per phenotype
  /// @param v_h2_gen_ Narrow-sense heritability per phenotype
  /// @param v_h2_env_ Environmental variance proportion per phenotype
  /// @param v_h2_vert_ Vertical transmission proportion per phenotype
  /// @param gen_cor_ Genetic correlation matrix
  /// @param env_cor_ Environmental correlation matrix
  /// @param v_rvert_pat_ Paternal / maternal vertical transmission ratio
  /// @param v_rvert_env_ Environmental / parental vertical transmission ratio
  ///
  /// @return Reference to this for method chaining
  SimulationConfig& phenome(
      std::size_t n_pheno_,
      const std::vector<std::string>& v_name_,
      std::vector<std::size_t> v_n_loc_,
      std::vector<double> v_h2_gen_,
      std::vector<double> v_h2_env_,
      std::vector<double> v_h2_vert_,
      std::optional<std::vector<double>> gen_cor_ = std::nullopt,
      std::optional<std::vector<double>> env_cor_ = std::nullopt,
      std::optional<std::vector<double>> v_rvert_pat_ = std::nullopt,
      std::optional<std::vector<double>> v_rvert_env_ = std::nullopt);

  /// @brief Enable random mating
  /// @return Reference to this for method chaining
  SimulationConfig& random_mating();

  /// @brief Enable assortative mating
  ///
  /// @param mate_cor_ Target mate phenotypic correlations
  /// @param tol_inf_ Tolerance for infinite values (default: 1e-7)
  /// @param n_itr_ Annealing iterations (default: 2e6)
  /// @param temp_init_ Initial temperature (default: 1.0)
  /// @param temp_decay_ Temperature decay (default: 0.999)
  /// @return Reference to this for method chaining
  SimulationConfig& assortative_mating(
      std::vector<double> mate_cor_,
      std::optional<double> tol_inf_ = std::nullopt,
      std::optional<std::size_t> n_itr_ = std::nullopt,
      std::optional<double> temp_init_ = std::nullopt,
      std::optional<double> temp_decay_ = std::nullopt);

  /// @brief Set metrics to compute each generation
  ///
  /// @param metric_specs_ Metric specifications
  /// @return Reference to this for method chaining
  SimulationConfig& metrics(std::vector<MetricSpec> metric_specs_);

  std::size_t n_gen;              ///< Number of generations
  std::size_t n_ind;              ///< Population size
  std::filesystem::path out_dir;  ///< Output directory
  std::uint64_t rng_seed;         ///< RNG seed

  std::size_t n_loc;          ///< Number of loci
  std::vector<double> v_maf;  ///< Minor allele frequencies
  std::vector<double> v_rec;  ///< Recombination probabilities
  std::vector<double> v_mut;  ///< Mutation probabilities

  std::size_t n_pheno;               ///< Number of phenotypes
  std::vector<std::string> v_name;   ///< Phenotype names
  std::vector<std::size_t> v_n_loc;  ///< Loci per phenotype
  std::vector<double> v_h2_gen;      ///< Narrow-sense heritabilities
  std::vector<double> v_h2_env;      ///< Environmental variance proportions
  std::vector<double> v_h2_vert;     ///< Vertical transmission proportions
  std::vector<double> gen_cor;       ///< Genetic correlation matrix
  std::vector<double> env_cor;       ///< Environmental correlation matrix
  std::vector<double> v_rvert_pat;   ///< Paternal transmission ratios
  std::vector<double> v_rvert_env;   ///< Nurture / environment ratios

  MatingType mating_type;        ///< Mating model type
  std::vector<double> mate_cor;  ///< Target mate correlations
  double tol_inf = 1e-7;         ///< Tolerance for infinite values
  std::size_t n_itr = 2e6;       ///< Annealing iterations
  double temp_init = 1.0;        ///< Initial annealing temperature
  double temp_decay = 0.99;      ///< Temperature decay factor

  std::vector<MetricSpec> specs;  ///< Metric specifications
  bool require_lat;               ///< Whether latent phenotypes are required
};

}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_CONFIG_H
