#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>

#include <Eigen/Dense>

namespace amsim {

class AssortativeMating {
 public:
  AssortativeMating(const AssortativeMating&) = default;
  AssortativeMating(AssortativeMating&&) = default;
  AssortativeMating& operator=(const AssortativeMating&) = default;
  AssortativeMating& operator=(AssortativeMating&&) = default;
  explicit AssortativeMating(const Params& params)
      : n_sex_(params.geno.n_ind / 2),
        n_pheno_(params.pheno.n_pheno),
        mate_cor_(std::move(params.mate.mate_cor)),
        match_cur_(n_sex_),
        match_opt_(n_sex_),
        max_itr_(params.mate.max_itr),
        temp_init_(params.mate.temp_init),
        temp_decay_(params.mate.temp_decay),
        tol_inf_(params.mate.tol_inf) {
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

}
