#include <amsim/state.h>
#include <amsim/transform.h>

namespace amsim::phenome {

void ScorePhenotypes::scoreGenetic(State& state) {
  if (state.geno().view() != genome::HaploView::LocusMajor)
    throw std::runtime_error("phenotype scoring requires loc-major view");

  auto& geno_buf = state.geno();
  auto pheno_gen_buf = state.pheno()(Component::Genetic);

  for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno) {
    const std::vector<std::size_t>& pheno_loc = pheno_loc_[pheno];
    const Eigen::VectorXd& pheno_effects = pheno_effects_[pheno];
    const std::size_t n_causal = pheno_loc.size();

    for (std::size_t ind_tile = 0; ind_tile < n_ind_; ind_tile += IND_TILE) {
      std::size_t tile_size = std::min(IND_TILE, n_ind_ - ind_tile);

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

void ScorePhenotypes::scoreEnvironmental(State& state) {
  auto pheno_env_buf = state.pheno()(Component::Environmental);
  rng::NormalPolar::fill(
      state.pheno()(Component::Environmental).data(), n_ind_ * n_pheno_);
  pheno_env_buf = pheno_env_buf * env_chol_.transpose();

  for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
    pheno_env_buf.col(pheno).noalias() =
        std::sqrt(h2_env_[pheno]) * pheno_env_buf.col(pheno);
}

void ScorePhenotypes::scoreNurture(State& state) {
  // this is the basic gist, of course we need to add scaling and both parents
  // for each child, but that isn't too hard since we have the matching :-D
  //
  // hard to compute everything at once, but we can do each phenotype at a time
  if (state.gen == 0) {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno) {
      auto pheno_nur_buf = state.pheno()(pheno, Component::Nurture);
      rng::NormalPolar::fill(pheno_nur_buf.data(), n_ind_);
      pheno_nur_buf *= std::sqrt(h2_nur_[pheno]);
    }
  } else {
    auto pheno_gen_buf = state.pheno()(Component::Genetic);
    auto pheno_par_gen_buf = state.pheno_par()(Component::Genetic);
    state.pheno()(Component::Nurture) = pheno_par_gen_buf - pheno_gen_buf;
  }
}

void ScorePhenotypes::scoreTotal(State& state) {
  state.pheno()(Component::Total).noalias() =
      state.pheno()(Component::Genetic) +
      state.pheno()(Component::Environmental);
}

void ScorePhenotypes::operator()(State& state) {
  scoreGenetic(state);
  scoreEnvironmental(state);
  scoreTotal(state);
}

}  // namespace amsim::phenome
