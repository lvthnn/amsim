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

class Simulation {
 public:
  explicit Simulation(
      const SimulationConfig& config,
      std::optional<std::filesystem::path> out_dir_ = std::nullopt,
      std::optional<std::uint64_t> rng_seed_ = std::nullopt);

  Simulation(Simulation&& other) noexcept
      : n_gen(other.n_gen),
        n_ind(other.n_ind),
        n_loc(other.n_loc),
        n_pheno(other.n_pheno),
        pheno_names(std::move(other.pheno_names)),
        out_dir(std::move(other.out_dir)),
        rng_(std::move(other.rng_)),
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

  const std::size_t n_gen;
  const std::size_t n_ind;
  const std::size_t n_loc;
  const std::size_t n_pheno;
  const std::vector<std::string> pheno_names;

  const std::filesystem::path out_dir;

  void run();

 private:
  rng::Xoshiro256ss rng_;
  Genome genome_;
  PhenoArch arch_;
  PhenoBuf buf_;
  PhenotypeList phenotypes_;
  AssortativeModel model_;
  SimulationContext ctx_;
  std::vector<Metric> metrics_;
  std::vector<std::unique_ptr<std::ofstream>> streams_;

  void stream(std::size_t gen);
};

void run_simulations(
    const SimulationConfig& config,
    std::size_t n_replicates,
    std::size_t n_threads,
    bool summarise = false,
    bool log_file = true,
    LogLevel log_level = LogLevel::INFO);
}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_H
