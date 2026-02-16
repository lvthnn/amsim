#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>

#include <ranges>
#include <vector>

namespace amsim {

class ScorePhenotypes {
 public:
  explicit ScorePhenotypes(const Params& params)
      : n_ind_(params.geno.n_ind),
        n_sex_(n_ind_ / 2),
        n_pheno_(params.pheno.n_pheno),
        pheno_effects_(params.pheno.pheno_effects),
        pheno_loc_(params.pheno.pheno_loc),
        h2_gen_(params.pheno.h2_gen),
        h2_env_(params.pheno.h2_env),
        h2_nur_(params.pheno.h2_nur),
        rnur_pat_(params.pheno.rnur_pat),
        rnur_env_(params.pheno.rnur_env),
        env_cor_(params.pheno.env_cor),
        nur_scale_(n_pheno_) {
    env_chol_ = env_cor_.llt().matrixL();

    std::size_t max_size = std::ranges::max(
        pheno_loc_ | std::views::transform(&std::vector<std::size_t>::size));

    gen_tile_.resize(IndTile, max_size);
  };

  void operator()(State& state);

 private:
  const std::size_t n_ind_;
  const std::size_t n_sex_;
  const std::size_t n_pheno_;

  const std::vector<Eigen::VectorXd>& pheno_effects_;
  const std::vector<std::vector<std::size_t>>& pheno_loc_;
  const Eigen::VectorXd& h2_gen_;
  const Eigen::VectorXd& h2_env_;
  const Eigen::VectorXd& h2_nur_;
  const Eigen::VectorXd& rnur_pat_;
  const Eigen::VectorXd& rnur_env_;

  static constexpr std::size_t IndTile = 512;
  Eigen::MatrixXd gen_tile_;
  Eigen::MatrixXd env_chol_;
  const Eigen::MatrixXd& env_cor_;
  Eigen::VectorXd nur_scale_;

  void scoreGenetic(State& state);
  void scoreEnvironmental(State& state);
  void scoreNurture(State& state);
  static void scoreTotal(State& state);
};

} // namespace amsim
