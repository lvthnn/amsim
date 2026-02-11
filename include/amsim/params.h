#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace amsim {

namespace genome {
/// ---------------------------------------------------------------------------
/// Genome parameters
/// ---------------------------------------------------------------------------
struct GenomeParams {
  std::size_t n_ind;
  std::size_t n_loc;
  Eigen::VectorXd v_maf;
  Eigen::VectorXd v_rec;
  Eigen::VectorXd v_mut;
};
}  // namespace genome

namespace phenome {
/// ---------------------------------------------------------------------------
/// Phenome parameters
/// ---------------------------------------------------------------------------
struct PhenomeParams {
  std::size_t n_pheno;              ///< Number of phenotypes
  std::vector<std::string> names;   ///< Vector of phenotype names
  std::vector<std::size_t> n_locs;  ///< Vector of number of loci per phenotype

  std::vector<Eigen::VectorXd> pheno_effects;  ///< Matrix of effect vectors
  std::vector<std::vector<std::size_t>> pheno_loc;  ///< Phenotype causal loci

  Eigen::VectorXd h2_gen;  ///< Vector of genetic component variances
  Eigen::VectorXd h2_env;  ///< Vector of environmental component variances
  Eigen::VectorXd h2_nur;  ///< Vector of nurture component variances

  Eigen::MatrixXd gen_cor;  ///< Genetic component correlation matrix
  Eigen::MatrixXd env_cor;  ///< Environmental component correlation matrix

  Eigen::VectorXd rnur_pat;  ///< Paternal ratio in nurture effect
  Eigen::VectorXd rnur_env;  ///< Environmental ratio in nurture effect
};
}  // namespace phenome

namespace mating {
/// ---------------------------------------------------------------------------
/// Mating parameters
/// ---------------------------------------------------------------------------
struct MatingParams {
  Eigen::MatrixXd mate_cor;
  double tol_inf;
  std::size_t max_itr;
  double temp_init;
  double temp_decay;
};
}  // namespace mating

/// ---------------------------------------------------------------------------
/// Simulation parameters
/// ---------------------------------------------------------------------------
struct SimulationParams {
  std::size_t n_reps;
  std::size_t n_threads;
  std::uint64_t rng_seed;
  std::filesystem::path out_dir;
};

// Lightweight object to hold all parameters
struct Params {
  genome::GenomeParams geno;
  phenome::PhenomeParams pheno;
  mating::MatingParams mate;
  SimulationParams sim;
};
}  // namespace amsim
