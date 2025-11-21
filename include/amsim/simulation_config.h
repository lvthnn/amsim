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

namespace amsim {

struct SimulationConfig {
  SimulationConfig() = default;

  SimulationConfig& simulation(
      std::size_t n_gen_,
      std::size_t n_ind_,
      const std::string& out_dir_,
      std::optional<std::uint64_t> rng_seed_);

  SimulationConfig& genome(
      std::size_t n_loc_,
      std::vector<double> v_maf_,
      std::vector<double> v_rec_,
      std::vector<double> v_mut_);

  SimulationConfig& phenome(
      std::size_t n_pheno_,
      const std::vector<std::string>& v_name_,
      std::vector<std::size_t> v_n_loc_,
      std::vector<double> v_h2_gen_,
      std::vector<double> v_h2_env_,
      std::vector<double> v_h2_vert_,
      std::vector<double> gen_cor_,
      std::vector<double> env_cor_);

  SimulationConfig& random_mating();

  SimulationConfig& assortative_mating(
      std::vector<double> mate_cor_,
      std::optional<double> tol_inf_ = 1e-7,
      std::optional<std::size_t> n_itr_ = 2e6,
      std::optional<double> temp_init_ = 1.0,
      std::optional<double> temp_decay_ = 0.999);

  SimulationConfig& metrics(std::vector<MetricSpec> metric_specs_);

  std::size_t n_gen;
  std::size_t n_ind;
  std::filesystem::path out_dir;
  std::uint64_t rng_seed;

  std::size_t n_loc;
  std::vector<double> v_maf;
  std::vector<double> v_rec;
  std::vector<double> v_mut;

  std::size_t n_pheno;
  std::vector<std::string> v_name;
  std::vector<std::size_t> v_n_loc;
  std::vector<double> v_h2_gen;
  std::vector<double> v_h2_env;
  std::vector<double> v_h2_vert;
  std::vector<double> gen_cor;
  std::vector<double> env_cor;

  MatingType mating_type;
  std::vector<double> mate_cor;
  double tol_inf;
  std::size_t n_itr;
  double temp_init;
  double temp_decay;

  std::vector<MetricSpec> specs;
  bool require_lat;
};

}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_CONFIG_H
