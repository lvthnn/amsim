#ifndef AMSIMCPP_SIMULATION_H
#define AMSIMCPP_SIMULATION_H

#include <amsim/genome.h>
#include <amsim/log_level.h>
#include <amsim/mating.h>
#include <amsim/metricspec.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/simulation_config.h>
#include <amsim/simulation_results.h>

#include <filesystem>
#include <fstream>
#include <vector>

namespace amsim {

/// @brief Forward-time assortative mating simulation
///
/// Simulation owns the genome, phenotypes, mating model, and metrics,
/// advancing the population forward through time.
class Simulation {
 public:
  /// @brief Construct a Simulation from configuration
  ///
  /// @param config Simulation configuration
  /// @param out_dir_ Optional output directory override
  /// @param rng_seed_ Optional RNG seed override
  explicit Simulation(
      const SimulationConfig& config,
      std::optional<std::filesystem::path> out_dir_ = std::nullopt,
      std::optional<std::uint64_t> rng_seed_ = std::nullopt);

  /// @brief Move constructor
  Simulation(Simulation&& other) noexcept
      : n_gen(other.n_gen),
        n_ind(other.n_ind),
        n_loc(other.n_loc),
        n_pheno(other.n_pheno),
        pheno_names(std::move(other.pheno_names)),
        out_dir(std::move(other.out_dir)),
        genome_(std::move(other.genome_)),
        arch_(std::move(other.arch_)),
        buf_(std::move(other.buf_)),
        phenotypes_(std::move(other.phenotypes_)),
        model_(std::move(other.model_)),
        ctx_(genome_, arch_, buf_, phenotypes_, model_),
        metrics_(std::move(other.metrics_)),
        streams_(std::move(other.streams_)) {}

  Simulation(const Simulation&) = delete;
  Simulation& operator=(const Simulation&) = delete;
  Simulation& operator=(Simulation&&) = delete;

  const std::size_t n_gen;                     ///< Number of generations
  const std::size_t n_ind;                     ///< Population size
  const std::size_t n_loc;                     ///< Number of loci
  const std::size_t n_pheno;                   ///< Number of phenotypes
  const std::vector<std::string> pheno_names;  ///< Phenotype names

  const std::filesystem::path out_dir;  ///< Output directory

  /// @brief Run the simulation
  void run();

 private:
  Genome genome_;                ///< Genome state
  PhenoArch arch_;               ///< Phenotype architecture
  PhenoBuf buf_;                 ///< Phenotype buffer
  PhenotypeList phenotypes_;     ///< Phenotype instances
  AssortativeModel model_;       ///< Mating model
  SimulationContext ctx_;        ///< Simulation context
  std::vector<Metric> metrics_;  ///< Metrics to compute
  std::vector<std::unique_ptr<std::ofstream>> streams_;  ///< Output streams

  /// @brief Write metrics for current generation to output streams
  /// @param gen Generation number
  void stream(std::size_t gen);
};

/// @brief Shuffle a seed based on the replicate simulation ID
///
/// @param rng_seed The base RNG seed used by the simulation
/// @param rep_id The replication ID used to shuffle the base seed
///
/// @return The shuffled seed.
std::uint64_t shuffle_seed(std::uint64_t rng_seed, std::size_t rep_id);

/// @brief Run a single simulation
///
/// @param config Simulation configuration
/// @param out_dir_ Overwrite configuration output directory
/// @param rng_seed_ Overwrite configuration random seed
/// @param log_file Whether to write log to file
/// @param log_level Logging verbosity level
void run_simulation(
    const SimulationConfig& config,
    std::optional<std::filesystem::path> out_dir_,
    bool log_file = false,
    LogLevel log_level = LogLevel::INFO);

/// @brief Run multiple simulation replicates in parallel
///
/// @param config Simulation configuration
/// @param n_replicates Number of replicate simulations
/// @param n_threads Number of parallel threads
/// @param summarise Whether to summarise results across replicates
/// @param log_file Whether to write log to file
/// @param log_level Logging verbosity level
void run_simulations(
    const SimulationConfig& config,
    std::size_t n_replicates,
    std::size_t n_threads,
    bool summarise = false,
    bool log_file = true,
    LogLevel log_level = LogLevel::INFO);
}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_H
