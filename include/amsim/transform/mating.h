// This file is part of amsim, copyright (C) 2025-2026 Kári Hlynsson.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace amsim {

class RandomMating {
 public:
  explicit RandomMating(const Params& params)
    : n_sex_(params.global.n_ind / 2),
      match_cur_(n_sex_) {}

  void operator()(State& state);

 private:
  std::size_t n_sex_;
  Matching match_cur_;

  void randomiseMatching();
};

inline void RandomMating::randomiseMatching() {
  std::iota(match_cur_.begin(), match_cur_.end(), 0);
  for (std::size_t el = 0; el < n_sex_; ++el) {
    std::size_t le = rng::UniformIntRange::sample(el, n_sex_);
    std::swap(match_cur_[el], match_cur_[le]);
  }
}

inline void RandomMating::operator()(State& state) {
  randomiseMatching();
  state.matching() = match_cur_;
  for (std::size_t ind = 0; ind < n_sex_; ++ind)
    state.inv_matching()[state.matching()[ind]] = ind;
}

class AssortativeMating {
 public:
  AssortativeMating(const AssortativeMating&) = default;
  AssortativeMating(AssortativeMating&&) = default;
  AssortativeMating& operator=(const AssortativeMating&) = default;
  AssortativeMating& operator=(AssortativeMating&&) = default;
  explicit AssortativeMating(const Params& params)
      : n_sex_(params.global.n_ind / 2),
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

inline void AssortativeMating::randomiseMatching() {
  std::iota(match_cur_.begin(), match_cur_.end(), 0);
  for (std::size_t el = 0; el < n_sex_; ++el) {
    std::size_t le = rng::UniformIntRange::sample(el, n_sex_);
    std::swap(match_cur_[el], match_cur_[le]);
  }
}

inline void AssortativeMating::latentMatching() {
  // compute rank-one phenotype projections
  Eigen::VectorXd latent_male = std_male_ * mate_cor_U1_;
  Eigen::VectorXd latent_female = std_female_ * mate_cor_V1_;
  Eigen::VectorXd fuzz(n_sex_);

  // add normal error to order-match with correlation mate_cor_S1_
  rng::NormalPolar::fill(fuzz.data(), n_sex_);
  double latent_noise = std::sqrt(1 - (mate_cor_S1_ * mate_cor_S1_));
  latent_female.noalias() += latent_noise * fuzz;

  // order mates, with females ranked on fuzzed projection score
  std::vector<std::size_t> idx_male = utils::order(latent_male);
  std::vector<std::size_t> idx_female = utils::order(latent_female);

  // assign mates according to order in rank-one projection
  for (std::size_t rank = 0; rank < n_sex_; ++rank)
    match_cur_[idx_male[rank]] = idx_female[rank];
}

inline void AssortativeMating::computeCor() {
  cor_ = Eigen::MatrixXd::Zero(n_pheno_, n_pheno_);
  for (std::size_t pair = 0; pair < n_sex_; ++pair)
    cor_.noalias() +=
        std_male_.row(pair).transpose() * std_female_.row(match_cur_[pair]);
  cor_ /= static_cast<double>(n_sex_);
}

inline void AssortativeMating::proposeState() {
  i0_ = amsim::rng::UniformIntRange::sample(n_sex_);
  i1_ = amsim::rng::UniformIntRange::sample(n_sex_);
  while (i0_ == i1_) i1_ = amsim::rng::UniformIntRange::sample(n_sex_);

  // compute the signed update matrix
  delta_cur_ =
      (std_male_.row(i0_) - std_male_.row(i1_)).transpose() *
      (std_female_.row(match_cur_[i1_]) - std_female_.row(match_cur_[i0_]));
  delta_cur_ /= n_sex_;

  // compute the energy differential
  alpha_cur_ =
      delta_cur_.squaredNorm() + (2 * delta_cur_.cwiseProduct(ell_cur_).sum());
}

inline void AssortativeMating::updateState() {
  // compute acceptance probability
  double acc_prob = std::min(1.0, std::exp(-alpha_cur_ / temp_cur_));
  double u = rng::UniformRange::sample(1.0);

  // accept
  if (u < acc_prob) {
    std::swap(match_cur_[i0_], match_cur_[i1_]);
    ell_cur_ += delta_cur_;

    double err_l2_cur = ell_cur_.norm();
    if (err_l2_cur < err_l2_opt_) {
      ell_opt_ = ell_cur_;
      err_l2_opt_ = err_l2_cur;
      err_linfty_opt_ = ell_cur_.array().abs().maxCoeff();
      match_opt_ = match_cur_;
    }
  }
}

inline void AssortativeMating::operator()(State& state) {
  const Eigen::MatrixXd& pheno_male = state.pheno().male();
  const Eigen::MatrixXd& pheno_female = state.pheno().female();

  // generate a random matching using Fisher-Yates shuffling
  randomiseMatching();

  // if max_itr_ has been assigned zero, we are doing random mating
  if (max_itr_ == 0) {
    n_itr_ = 0;
    return;
  }

  n_itr_ = max_itr_;

  // standardise the supplied matrices
  std_male_ = utils::standardise(pheno_male);
  std_female_ = utils::standardise(pheno_female);

  latentMatching();
  computeCor();

  // initialise values for this run
  match_opt_ = match_cur_;
  temp_cur_ = temp_init_;
  ell_cur_ = cor_ - mate_cor_;
  ell_opt_ = ell_cur_;
  err_l2_opt_ = ell_cur_.norm();
  err_linfty_opt_ = ell_cur_.array().abs().maxCoeff();

  // annealing routine
  for (std::size_t itr = 0; itr < max_itr_; ++itr) {
    proposeState();
    updateState();
    temp_cur_ *= temp_decay_;

    if (err_linfty_opt_ < tol_inf_) {
      n_itr_ = itr;
      break;
    }
  }

  n_itr_ = max_itr_;
  cor_ = ell_opt_ + mate_cor_;
  state.matching() = match_opt_;

  for (std::size_t ind = 0; ind < n_sex_; ++ind)
    state.inv_matching()[state.matching()[ind]] = ind;
}

} // namespace amsim
