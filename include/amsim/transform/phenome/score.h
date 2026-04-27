#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>

#include <ranges>
#include <vector>
#include <stdexcept>

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
        h2_vert_(params.pheno.h2_vert),
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
  const Eigen::VectorXd& h2_vert_;
  const Eigen::VectorXd& rnur_pat_;
  const Eigen::VectorXd& rnur_env_;

  static constexpr std::size_t IndTile = 512;
  Eigen::MatrixXd gen_tile_;
  Eigen::MatrixXd env_chol_;
  const Eigen::MatrixXd& env_cor_;
  Eigen::VectorXd nur_scale_;

  void scoreGenetic(State& state);
  void scoreEnvironmental(State& state);
  void scoreVertical(State& state);
  static void scoreTotal(State& state);
};

inline void ScorePhenotypes::scoreGenetic(State& state) {
  if (state.geno().view() != HaploView::LocusMajor)
    throw std::runtime_error("phenotype scoring requires loc-major view");

  auto& geno_buf = state.geno();
  auto pheno_gen_buf = state.pheno()(Component::Genetic);

  for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno) {
    const std::vector<std::size_t>& pheno_loc = pheno_loc_[pheno];
    const Eigen::VectorXd& pheno_effects = pheno_effects_[pheno];
    const std::size_t n_causal = pheno_loc.size();

    for (std::size_t ind_tile = 0; ind_tile < n_ind_; ind_tile += IndTile) {
      std::size_t tile_size = std::min(IndTile, n_ind_ - ind_tile);

      // decompress the standardised genotype matrix in manageable tiles
      geno_buf.decompress(
          ind_tile, ind_tile + tile_size, pheno_loc, gen_tile_, true);

      // update the genetic component buffer in-place
      pheno_gen_buf.col(pheno).segment(ind_tile, tile_size).noalias() =
          std::sqrt(h2_gen_[pheno]) *
          gen_tile_.topRows(tile_size).leftCols(n_causal) * pheno_effects;
    }
  }
}

inline void ScorePhenotypes::scoreEnvironmental(State& state) {
  auto pheno_env_buf = state.pheno()(Component::Environmental);
  rng::NormalPolar::fill(
      state.pheno()(Component::Environmental).data(), n_ind_ * n_pheno_);
  pheno_env_buf = pheno_env_buf * env_chol_.transpose();

  for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
    pheno_env_buf.col(pheno).noalias() =
        std::sqrt(h2_env_[pheno]) * pheno_env_buf.col(pheno);
}

inline void ScorePhenotypes::scoreVertical(State& state) {
  if (h2_vert_.isZero(0)) return;
}

inline void ScorePhenotypes::scoreTotal(State& state) {
  state.pheno()(Component::Total).noalias() =
      state.pheno()(Component::Genetic) +
      state.pheno()(Component::Environmental) +
      state.pheno()(Component::Vertical);
}

inline void ScorePhenotypes::operator()(State& state) {
  scoreGenetic(state);
  scoreEnvironmental(state);
  scoreNurture(state);
  scoreTotal(state);
}

}  // namespace amsim
