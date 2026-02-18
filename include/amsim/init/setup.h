#pragma once

#include <amsim/core.h>
#include <amsim/estimate.h>
#include <amsim/sample.h>

#include <Eigen/Dense>
#include <cstddef>
#include <string>
#include <vector>
namespace amsim {

namespace details {

inline Eigen::VectorXd expand(
    std::variant<double, Eigen::VectorXd> val, std::size_t size) {
  if (std::holds_alternative<double>(val))
    return Eigen::VectorXd::Constant(size, std::get<double>(val));

  Eigen::VectorXd res = std::get<Eigen::VectorXd>(val);

  if (static_cast<std::size_t>(res.size()) != size)
    throw std::runtime_error(
        "incorrect vector size; got " + std::to_string(res.size()) +
        ", expected " + std::to_string(size));

  return res;
}

}  // namespace details

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

  std::vector<PopulationEstimator> estimators;  // population-wide estimators
  std::vector<SampleSpec> samples;              // subpopulation estimators

  std::size_t n_generations;
  std::filesystem::path output_dir;
  std::optional<std::uint64_t> random_seed;

  LogLevel log_level;
  bool log_file = false;
};

inline PhenomeParams build_pheno_params(const Simulation& simulation) {
  std::size_t n_pheno = simulation.phenotypes.size();

  std::vector<std::string> names(n_pheno);
  std::vector<std::size_t> n_locs(n_pheno);
  std::unordered_map<std::string, std::size_t> pheno_ids(n_pheno);

  std::vector<Eigen::VectorXd> pheno_effects(n_pheno);
  std::vector<std::vector<std::size_t>> pheno_loc(n_pheno);

  Eigen::VectorXd h2_gen(n_pheno);
  Eigen::VectorXd h2_env(n_pheno);
  Eigen::VectorXd h2_nur(n_pheno);

  Eigen::VectorXd rnur_pat(n_pheno);
  Eigen::VectorXd rnur_env(n_pheno);

  for (std::size_t pheno = 0; pheno < n_pheno; ++pheno) {
    Phenotype pheno_data = simulation.phenotypes[pheno];

    names[pheno] = pheno_data.name;
    n_locs[pheno] = pheno_data.n_causal_loci;
    pheno_ids[pheno_data.name] = pheno;

    if (!pheno_data.effects.has_value())
      pheno_effects[pheno] = Eigen::VectorXd::Constant(
          pheno_data.n_causal_loci, 1.0 / std::sqrt(pheno_data.n_causal_loci));
    else
      pheno_effects[pheno] =
          pheno_data.effects->array() / pheno_data.effects->norm();

    if (pheno_data.causal_loci.has_value() &&
        simulation.genetic_component_cor.has_value())
      throw std::runtime_error(
          "specify either exclusively phenotype causal loci or genetic "
          "component correlation");

    if (!pheno_data.causal_loci.has_value() &&
        !simulation.genetic_component_cor.has_value())
      throw std::runtime_error(
          "require either causal loci or genetic component correlation to "
          "intialise phenotypes");

    if (pheno_data.causal_loci.has_value()) {
      if (pheno_data.causal_loci.value().size() != pheno_data.n_causal_loci)
        throw std::runtime_error(
            "phenotype " + pheno_data.name + "expected " +
            std::to_string(pheno_data.n_causal_loci) + " causal loci, got " +
            std::to_string(pheno_data.causal_loci.value().size()));
      pheno_loc[pheno] = pheno_data.causal_loci.value();
    } else
      pheno_loc[pheno].resize(pheno_data.n_causal_loci);

    h2_gen(pheno) = pheno_data.h2_genetic;
    h2_env(pheno) = pheno_data.h2_environmental;
    h2_nur(pheno) = pheno_data.h2_nurture;
    rnur_pat(pheno) = pheno_data.nurture_paternal_ratio;
    rnur_env(pheno) = pheno_data.nurture_environmental_ratio;
  }

  Eigen::MatrixXd gen_cor = simulation.genetic_component_cor.has_value()
                                ? simulation.genetic_component_cor.value()
                                : Eigen::MatrixXd::Identity(n_pheno, n_pheno);

  Eigen::MatrixXd env_cor = simulation.environmental_component_cor.has_value()
                                ? simulation.environmental_component_cor.value()
                                : Eigen::MatrixXd::Identity(n_pheno, n_pheno);

  return PhenomeParams{
      .n_pheno = n_pheno,
      .names = std::move(names),
      .n_locs = std::move(n_locs),
      .pheno_ids = std::move(pheno_ids),
      .pheno_effects = std::move(pheno_effects),
      .pheno_loc = std::move(pheno_loc),
      .h2_gen = h2_gen,
      .h2_env = h2_env,
      .h2_nur = h2_nur,
      .gen_cor = gen_cor,
      .env_cor = env_cor,
      .rnur_pat = rnur_pat,
      .rnur_env = rnur_env};
}

inline Params build_params(const Simulation& simulation) {
  GenomeParams geno = GenomeParams{
      .n_ind = simulation.n_individuals,
      .n_loc = simulation.genome.n_loci,
      .v_rec =
          details::expand(simulation.genome.v_rec, simulation.genome.n_loci),
      .v_maf =
          details::expand(simulation.genome.v_maf, simulation.genome.n_loci),
      .v_mut =
          details::expand(simulation.genome.v_mut, simulation.genome.n_loci)};

  PhenomeParams pheno = build_pheno_params(simulation);

  MatingParams mate = MatingParams{
      .mate_cor = std::move(simulation.mating.mate_cor),
      .tol_inf = simulation.mating.tolerance,
      .max_itr = simulation.mating.max_iterations,
      .temp_init = simulation.mating.initial_temperature,
      .temp_decay = simulation.mating.temperature_decay};

  SimulationParams sim = SimulationParams{
      .n_gens = simulation.n_generations,
      .out_dir = simulation.output_dir,
      .rng_seed = rng::auto_seed(simulation.random_seed)};

  return Params{
      .geno = std::move(geno),
      .pheno = std::move(pheno),
      .mate = std::move(mate),
      .sim = std::move(sim)};
}

}  // namespace amsim
