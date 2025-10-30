#ifndef AMSIMCPP_SIMULATION_BUILDER_H
#define AMSIMCPP_SIMULATION_BUILDER_H

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

namespace amsim {

class SimulationBuilder {
 public:
  SimulationBuilder();

  SimulationBuilder& simulation(
      std::size_t n_gen,
      std::size_t n_ind,
      std::string out_dir,
      std::uint64_t rng_seed,
      std::optional<std::uint64_t> rep_id = std::nullopt);

  SimulationBuilder& genome(
      std::size_t n_loc,
      std::vector<double> v_maf,
      std::vector<double> v_rec,
      std::vector<double> v_mut);

  SimulationBuilder& phenome(
      std::size_t n_pheno,
      std::vector<std::string> v_name,
      std::vector<std::size_t> v_n_loc,
      std::vector<double> v_h2_gen,
      std::vector<double> v_h2_env,
      std::vector<double> v_h2_vert,
      std::vector<double> gen_cor,
      std::vector<double> env_cor);

  SimulationBuilder& mating(
      MatingType type,
      std::optional<std::size_t> n_itr,
      std::optional<double> tmp_init,
      std::optional<double> tmp_decay,
      std::optional<std::vector<double>> mate_cor);

  SimulationBuilder& metrics(std::vector<MetricSpec> metric_specs);

  Simulation build();

 private:
  //--- OTHER PARAMETERS ----------------------------------------------------//
  SimulationStatus status_;
  std::optional<std::size_t> rep_id_;

  //--- SIMULATION PARAMETERS -----------------------------------------------//
  std::size_t n_gen_;
  std::size_t n_ind_;
  std::filesystem::path out_dir_;
  std::uint64_t rng_seed_;

  //--- GENOME PARAMETERS ---------------------------------------------------//
  std::size_t n_loc_;
  std::vector<double> v_maf_;
  std::vector<double> v_rec_;
  std::vector<double> v_mut_;

  //--- PHENOTYPE PARAMETERS ------------------------------------------------//
  std::size_t n_pheno_;
  std::vector<std::string> v_name_;
  std::vector<std::size_t> v_n_loc_;
  std::vector<double> v_h2_gen_;
  std::vector<double> v_h2_env_;
  std::vector<double> v_h2_vert_;
  std::vector<double> gen_cor_;
  std::vector<double> env_cor_;

  //--- MATING MODEL PARAMETERS ---------------------------------------------//
  std::size_t n_itr_;
  double tmp_init_;
  double tmp_decay_;
  std::vector<double> mate_cor_;

  //--- METRIC PARAMETERS ---------------------------------------------------//
  std::vector<MetricSpec> specs_;
  bool require_lat_ = false;
};

}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_BUILDER_H
