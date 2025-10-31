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
#include <amsim/simulation.h>
#include <amsim/simulation_status.h>

#include <filesystem>

namespace amsim {

struct SimulationConfig {
  SimulationConfig();

  SimulationConfig& simulation(
      std::size_t n_gen_,
      std::size_t n_ind_,
      std::string out_dir_,
      std::uint64_t rng_seed_);

  SimulationConfig& genome(
      std::size_t n_loc_,
      std::vector<double> v_maf_,
      std::vector<double> v_rec_,
      std::vector<double> v_mut_);

  SimulationConfig& phenome(
      std::size_t n_pheno_,
      std::vector<std::string> v_name_,
      std::vector<std::size_t> v_n_loc_,
      std::vector<double> v_h2_gen_,
      std::vector<double> v_h2_env_,
      std::vector<double> v_h2_vert_,
      std::vector<double> gen_cor_,
      std::vector<double> env_cor_);

  SimulationConfig& mating(
      MatingType type_,
      std::optional<std::size_t> n_itr_ = std::nullopt,
      std::optional<double> temp_init_ = std::nullopt,
      std::optional<double> temp_decay_ = std::nullopt,
      std::optional<std::vector<double>> mate_cor_ = std::nullopt);

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

  MatingType type;
  std::optional<std::size_t> n_itr;
  std::optional<double> temp_init;
  std::optional<double> temp_decay;
  std::optional<std::vector<double>> mate_cor;

  std::vector<MetricSpec> specs;
  bool require_lat;
};

}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_CONFIG_H
