#pragma once

#include <amsim/output/log_level.h>

#include <Eigen/Dense>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace amsim {

namespace genome {
/// ---------------------------------------------------------------------------
/// Genome parameters
/// ---------------------------------------------------------------------------
struct GenomeParams {
  GenomeParams(
      std::size_t n_ind_,
      std::size_t n_loc_,
      Eigen::VectorXd v_maf_,
      Eigen::VectorXd v_rec_,
      Eigen::VectorXd v_mut_,
      std::optional<Eigen::MatrixXd> ld_cor_ = std::nullopt,
      std::optional<std::size_t> ld_size_ = std::nullopt)
      : n_ind(n_ind_),
        n_loc(n_loc_),
        v_maf(std::move(v_maf_)),
        v_rec(std::move(v_rec_)),
        v_mut(std::move(v_mut_)),
        ld_cor(std::move(ld_cor_)),
        ld_size(std::move(ld_size_)) {};

  std::size_t n_ind;
  std::size_t n_loc;
  Eigen::VectorXd v_maf;
  Eigen::VectorXd v_rec;
  Eigen::VectorXd v_mut;

  std::optional<Eigen::MatrixXd> ld_cor;
  std::optional<std::size_t> ld_size;

  bool has_ld() const { return ld_cor.has_value(); }
};
} // namespace genome

namespace phenome {
/// ---------------------------------------------------------------------------
/// Phenome parameters
/// ---------------------------------------------------------------------------
struct PhenomeParams {
  explicit PhenomeParams(
      std::size_t n_pheno_,
      std::vector<std::string> names_,
      std::vector<std::size_t> n_locs_,
      std::optional<Eigen::VectorXd> h2_gen_ = std::nullopt,
      std::optional<Eigen::VectorXd> h2_env_ = std::nullopt,
      std::optional<Eigen::VectorXd> h2_nur_ = std::nullopt,
      std::optional<Eigen::MatrixXd> gen_cor_ = std::nullopt,
      std::optional<Eigen::MatrixXd> env_cor_ = std::nullopt,
      std::optional<Eigen::VectorXd> rnur_pat_ = std::nullopt,
      std::optional<Eigen::VectorXd> rnur_env_ = std::nullopt)
      : n_pheno(n_pheno_),
        names(std::move(names_)),
        n_locs(std::move(n_locs_)),
        h2_gen(
            h2_gen_ ? std::move(*h2_gen_)
                    : Eigen::VectorXd::Constant(n_pheno, 0.5)),
        h2_env(
            h2_env_ ? std::move(*h2_env_)
                    : Eigen::VectorXd::Constant(n_pheno, 0.5)),
        h2_nur(
            h2_nur_ ? std::move(*h2_nur_)
                    : Eigen::VectorXd::Constant(n_pheno, 0.0)),
        gen_cor(
            gen_cor_ ? std::move(*gen_cor_)
                     : Eigen::MatrixXd::Identity(n_pheno, n_pheno)),
        env_cor(
            env_cor_ ? std::move(*env_cor_)
                     : Eigen::MatrixXd::Identity(n_pheno, n_pheno)),
        rnur_pat(
            rnur_pat_ ? std::move(*rnur_pat_)
                      : Eigen::VectorXd::Constant(n_pheno, 0.5)),
        rnur_env(
            rnur_env_ ? std::move(*rnur_env_)
                      : Eigen::VectorXd::Constant(n_pheno, 0.5)) {};

  std::size_t n_pheno;             ///< Number of phenotypes
  std::vector<std::string> names;  ///< Vector of phenotype names
  std::vector<std::size_t> n_locs; ///< Vector of number of loci per phenotype

  Eigen::VectorXd h2_gen; ///< Vector of genetic component variances
  Eigen::VectorXd h2_env; ///< Vector of environmental component variances
  Eigen::VectorXd h2_nur; ///< Vector of nurture component variances

  Eigen::MatrixXd gen_cor; ///< Genetic component correlation matrix
  Eigen::MatrixXd env_cor; ///< Environmental component correlation matrix

  Eigen::MatrixXd pheno_effects;       ///< Matrix with effect column vectors 
  std::vector<std::size_t> causal_loc; ///< Vector of causal loci for phenotypes

  Eigen::VectorXd rnur_pat; ///< Paternal ratio in nurture effect
  Eigen::VectorXd rnur_env; ///< Environmental ratio in nurture effect
};
} // namespace phenome


namespace mating {
/// ---------------------------------------------------------------------------
/// Mating parameters
/// ---------------------------------------------------------------------------
struct MatingParams {
  explicit MatingParams(
      Eigen::MatrixXd mate_cor_,
      double tol_inf_ = 1e-7,
      std::size_t max_itr_ = 2e6,
      double temp_init_ = 1.0,
      double temp_decay_ = 0.99)
      : mate_cor(std::move(mate_cor_)),
        tol_inf(tol_inf_),
        max_itr(max_itr_),
        temp_init(temp_init_),
        temp_decay(temp_decay_) {};

  Eigen::MatrixXd mate_cor;
  Eigen::MatrixXd cor_U;
  Eigen::MatrixXd cor_S;
  Eigen::MatrixXd cor_VT;

  double tol_inf;
  std::size_t max_itr;
  double temp_init;
  double temp_decay;
};
} // namespace mating

/// ---------------------------------------------------------------------------
/// Simulation parameters
/// ---------------------------------------------------------------------------
struct SimulationParams {
  std::size_t n_reps;
  std::size_t n_threads;
  std::uint64_t rng_seed;

  std::filesystem::path out_dir;

  LogLevel log_level;
  bool log_to_file = false;
};

// Lightweight object to hold all parameters
struct Params {
  Params(
      genome::GenomeParams& geno_,
      phenome::PhenomeParams& pheno_,
      mating::MatingParams& mate_,
      SimulationParams& sim_)
      : geno(std::move(geno_)),
        pheno(std::move(pheno_)),
        mate(std::move(mate_)),
        sim(std::move(sim_)) {};

  genome::GenomeParams geno;
  phenome::PhenomeParams pheno;
  mating::MatingParams mate;
  SimulationParams sim;
};
}  // namespace amsim
