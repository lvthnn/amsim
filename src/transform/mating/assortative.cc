#include <amsim/core/rng.h>
#include <amsim/core/utils.h>
#include <amsim/transform/mating/assortative.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace amsim {

void AssortativeMating::randomiseMatching() {
  std::iota(match_cur_.begin(), match_cur_.end(), 0);
  for (std::size_t el = 0; el < n_sex_; ++el) {
    std::size_t le = rng::UniformIntRange::sample(el, n_sex_);
    std::swap(match_cur_[el], match_cur_[le]);
  }
}

void AssortativeMating::latentMatching() {
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

void AssortativeMating::computeCor() {
  cor_ = Eigen::MatrixXd::Zero(n_pheno_, n_pheno_);
  for (std::size_t pair = 0; pair < n_sex_; ++pair)
    cor_.noalias() +=
        std_male_.row(pair).transpose() * std_female_.row(match_cur_[pair]);
  cor_ /= static_cast<double>(n_sex_);
}

void AssortativeMating::proposeState() {
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
      delta_cur_.squaredNorm() + 2 * delta_cur_.cwiseProduct(ell_cur_).sum();
}

void AssortativeMating::updateState() {
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

void AssortativeMating::operator()(State& state) {
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
    state.inv_matching()[ind] = state.matching()[ind];
}

}  // namespace amsim
