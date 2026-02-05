#pragma once

#include <amsim/params.h>
#include <amsim/state.h>

#include <ranges>
#include <vector>

namespace amsim {

namespace genome {

// Class to update the genome
class UpdateGenome {
 public:
  explicit UpdateGenome(const Params& params)
      : n_ind_(params.geno.n_ind),
        n_sex_(n_ind_ / 2),
        v_rec_(params.geno.v_rec),
        v_mut_(params.geno.v_mut),
        bw_() {};
  void operator()(State& state);

 private:
  const std::size_t n_ind_;
  const std::size_t n_sex_;
  const Eigen::VectorXd& v_rec_;
  const Eigen::VectorXd& v_mut_;
  const double* ptr_rec_;
  const double* ptr_mut_;
  std::vector<std::uint64_t> transmit_chunk_;
  rng::BernoulliWord<16> bw_;

  std::array<std::uint64_t, 2> gamWord(std::uint64_t h0, std::uint64_t h1);

  void updateGenome(State& state);
  void updateNurture(State& state);
};

}  // namespace genome

namespace phenome {

// Scores phenotypes using GenoBuf and PhenoBuf
class ScorePhenotypes {
 public:
  explicit ScorePhenotypes(const Params& params)
      : n_ind_(params.geno.n_ind),
        n_pheno_(params.pheno.n_pheno),
        pheno_effects_(params.pheno.pheno_effects),
        pheno_loc_(params.pheno.causal_loc),
        h2_gen_(params.pheno.h2_gen),
        h2_env_(params.pheno.h2_env),
        h2_nur_(params.pheno.h2_nur),
        rnur_pat_(params.pheno.rnur_pat),
        rnur_env_(params.pheno.rnur_env),
        env_cor_(params.pheno.env_cor) {
    env_chol_ = env_cor_.llt().matrixL();

    std::size_t max_size = std::ranges::max(
        pheno_loc_ | std::views::transform(&std::vector<std::size_t>::size));

    gen_tile_.resize(IND_TILE, max_size);
  };

  void operator()(State& state);

 private:
  const std::size_t n_ind_;
  const std::size_t n_pheno_;

  const std::vector<Eigen::VectorXd>& pheno_effects_;
  const std::vector<std::vector<std::size_t>>& pheno_loc_;
  const Eigen::VectorXd& h2_gen_;
  const Eigen::VectorXd& h2_env_;
  const Eigen::VectorXd& h2_nur_;
  const Eigen::VectorXd& rnur_pat_;
  const Eigen::VectorXd& rnur_env_;

  static constexpr std::size_t IND_TILE = 512;
  Eigen::MatrixXd gen_tile_;
  Eigen::MatrixXd env_chol_;
  const Eigen::MatrixXd& env_cor_;

  void scoreGenetic(State& state);
  void scoreEnvironmental(State& state);
  void scoreTotal(State& state);
};

}  // namespace phenome

namespace mating {

class AssortativeMating {
 public:
  AssortativeMating(const AssortativeMating&) = default;
  AssortativeMating(AssortativeMating&&) = default;
  AssortativeMating& operator=(const AssortativeMating&) = default;
  AssortativeMating& operator=(AssortativeMating&&) = default;
  explicit AssortativeMating(const Params& p)
      : n_sex_(p.geno.n_ind / 2),
        n_pheno_(p.pheno.n_pheno),
        mate_cor_(std::move(p.mate.mate_cor)),
        match_cur_(n_sex_),
        match_opt_(n_sex_),
        max_itr_(p.mate.max_itr),
        temp_init_(p.mate.temp_init),
        temp_decay_(p.mate.temp_decay),
        tol_inf_(p.mate.tol_inf) {
    // initialise rank-one SVD components
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        mate_cor_, Eigen::ComputeThinU | Eigen::ComputeThinV);

    mate_cor_S1_ = svd.singularValues()(0);
    mate_cor_U1_ = svd.matrixU().col(0);
    mate_cor_V1_ = svd.matrixV().col(0);
  };

  void operator()(State& state);

  // optimal matching
  Matching matching() const { return match_opt_; };

  // best l2 error achieved by solver
  double l2() const { return err_l2_opt_; }

  // best linfty error achieved by solver
  double linfty() const { return err_linfty_opt_; }

  // best correlation structure achieved by solver
  Eigen::MatrixXd cor() const { return cor_; };

  // number of iterations performed by solver
  std::size_t n_itr() const { return n_itr_; };

 private:
  std::size_t n_sex_;
  std::size_t n_pheno_;
  Eigen::MatrixXd std_male_;
  Eigen::MatrixXd std_female_;
  Eigen::MatrixXd mate_cor_;

  double mate_cor_S1_;
  Eigen::VectorXd mate_cor_U1_;
  Eigen::VectorXd mate_cor_V1_;

  std::size_t i0_;
  std::size_t i1_;
  Matching match_cur_;
  Matching match_opt_;
  Eigen::MatrixXd delta_cur_;
  Eigen::MatrixXd ell_cur_;
  Eigen::MatrixXd ell_opt_;
  double alpha_cur_;
  double err_l2_opt_;
  double err_linfty_opt_;

  Eigen::MatrixXd cor_;
  std::size_t n_itr_;
  std::size_t max_itr_;

  double temp_init_;
  double temp_decay_;
  double temp_cur_;
  double tol_inf_;

  // shuffle the internal matching using Fisher-Yates
  void randomiseMatching();

  // match mates by rank-one latent score
  void latentMatching();

  // compute the cross-correlation matrix for the initial state
  void computeCor();

  // propose next state and compute acceptance probability
  void proposeState();

  // accept or reject the proposed step
  void updateState();
};

}  // namespace mating

}  // namespace amsim
