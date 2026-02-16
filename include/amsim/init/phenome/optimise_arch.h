#pragma once

#include <amsim/core/params.h>

#include <Eigen/Dense>
#include <vector>

namespace amsim {

// phenotype genetic component correlation optimiser
class OptimisePhenotypeArchitecture {
 public:
  explicit OptimisePhenotypeArchitecture(Params& params)
      : n_pheno_(params.pheno.n_pheno),
        n_loc_(params.geno.n_loc),
        n_locs_(params.pheno.n_locs),
        gen_cor_(params.pheno.gen_cor),
        pheno_effects_vec_(params.pheno.pheno_effects),
        pheno_loc_(params.pheno.pheno_loc),
        pheno_loc_complement_(n_pheno_),
        pheno_effects_(n_loc_, n_pheno_),
        diff_cur_(n_pheno_, n_pheno_) {}

  double err_frob() const { return err_frob_opt_; }

  Eigen::MatrixXd expected() const { return expected_; }

  void operator()();

 private:
  const std::size_t n_pheno_;
  const std::size_t n_loc_;
  const std::vector<std::size_t>& n_locs_;
  const Eigen::MatrixXd& gen_cor_;
  const std::vector<Eigen::VectorXd>& pheno_effects_vec_;
  std::vector<std::vector<std::size_t>>& pheno_loc_;
  std::vector<std::vector<std::size_t>> pheno_loc_complement_;
  Eigen::MatrixXd pheno_effects_;

  std::size_t max_itr_ = 2e6;
  double tol_l2_ = 1e-15;
  double temp_init_ = 1.0;
  double temp_decay_ = 0.9999;
  std::size_t n_itr_;

  Eigen::MatrixXd expected_;
  Eigen::MatrixXd diff_cur_;
  Eigen::VectorXd u_cur_;
  double d_cur_;
  double temp_cur_;

  std::size_t pheno_cur_;
  std::size_t loc_ind_cur_;
  std::size_t loc_ind_prop_;
  std::size_t loc_cur_;
  std::size_t loc_prop_;

  double delta_cur_;
  double err_frob_;
  double err_frob_opt_;

  // generate the initial state
  void randomState();

  // compute the initial objective
  void computeInitObjective();

  void printObjective();

  // compute the expected correlation
  void computeExpected();

  // generate a proposal state (phenotype + switch) and compute update
  void proposeState();

  // metropolis accept step and update
  void updateState();
};

}  // namespace amsim
