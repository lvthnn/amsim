#pragma once

#include <amsim/output/estimator.h>
#include <amsim/output/sample.h>
#include <amsim/params.h>

#include <Eigen/Dense>
#include <cstddef>
#include <string>
#include <vector>

namespace amsim {

struct Genome {
  std::size_t n_loci;

  std::variant<double, Eigen::VectorXd> v_rec = 0.5;
  std::variant<double, Eigen::VectorXd> v_maf = 0.5;
  std::variant<double, Eigen::VectorXd> v_mut = 0.0;
};

struct Phenotype {
  std::string name;
  std::size_t n_causal_loci;
  std::unordered_map<std::string, std::size_t> ids;

  std::optional<Eigen::VectorXd> effects;
  std::optional<std::vector<std::size_t>> causal_loci;

  double h2_genetic = 0.5;
  double h2_environmental = 0.5;
  double h2_nurture = 0.0;

  double nurture_paternal_ratio = 0.5;
  double nurture_environmental_ratio = 0.5;
};

struct Mating {
  Eigen::MatrixXd mate_cor;

  double tolerance = 1e-7;
  std::size_t max_iterations = 2000000;
  double initial_temperature = 1.0;
  double temperature_decay = 0.99;
};

struct Simulation {
  std::size_t n_individuals;

  Genome genome;
  std::vector<Phenotype> phenotypes;
  std::optional<Eigen::MatrixXd> genetic_component_cor;
  std::optional<Eigen::MatrixXd> environmental_component_cor;

  Mating mating;

  std::vector<PopulationEstimator> estimators; // population-wide estimators
  std::vector<SampleSpec> samples; // subpopulation estimators

  std::size_t n_replicates = 1;
  std::size_t n_threads = 1;
  std::filesystem::path output_dir;
  std::optional<std::uint64_t> random_seed;
};

Params build_params(const Simulation& simulation);

}  // namespace amsim
