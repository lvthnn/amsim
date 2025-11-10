#ifndef AMSIMCPP_SIMULATION_H
#define AMSIMCPP_SIMULATION_H

#pragma once

#include <filesystem>
#include <fstream>
#include <vector>

#include <amsim/genome.h>
#include <amsim/log_level.h>
#include <amsim/mating.h>
#include <amsim/metricspec.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/simulation_config.h>
#include <amsim/simulation_results.h>

namespace amsim {

class Simulation {
 public:
  Simulation(
      const SimulationConfig& config,
      std::optional<std::filesystem::path> out_dir_ = std::nullopt,
      std::optional<std::uint64_t> rng_seed_ = std::nullopt);

  Simulation(
      std::size_t n_gen,
      std::size_t n_ind,
      std::filesystem::path out_dir,
      std::uint64_t rng_seed,
      std::size_t n_loc,
      std::vector<double> v_maf,
      std::vector<double> v_rec,
      std::vector<double> v_mut,
      std::size_t n_pheno,
      std::vector<std::string> v_name,
      std::vector<std::size_t> v_n_loc,
      std::vector<double> v_h2_gen,
      std::vector<double> v_h2_env,
      std::vector<double> v_h2_vert,
      std::vector<double> gen_cor,
      std::vector<double> env_cor,
      std::vector<double> mate_cor,
      std::size_t n_itr,
      double tmp_init,
      double tmp_decay,
      std::vector<MetricSpec> specs,
      bool require_latent);

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

  void stream_(std::size_t gen);
};

void run_simulation(
    const SimulationConfig& config,
    std::optional<std::size_t> rep_id = std::nullopt);

void run_simulations(
    const SimulationConfig& config,
    std::size_t n_replicates,
    std::size_t n_threads,
    bool summarise = true,
    LogLevel level = LogLevel::INFO);

}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_H
