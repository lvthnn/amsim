#include <amsim/state.h>
#include <amsim/transform.h>

namespace amsim::phenome {

void ScorePhenotypes::scoreGenetic(State& state) {
  if (state.geno.view() != genome::HaploView::LOC_MAJOR)
    throw std::runtime_error("phenotype scoring requires loc-major view");

  auto& geno_buf = state.geno;
  auto pheno_gen_buf = state.pheno(ComponentType::GENETIC);

  for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno) {
    const std::vector<std::size_t>& pheno_loc = pheno_loc_[pheno];
    const Eigen::VectorXd& pheno_effects = pheno_effects_[pheno];
    const std::size_t n_causal = pheno_loc.size();

    for (std::size_t ind_tile = 0; ind_tile < n_ind_; ind_tile += IND_TILE) {
      std::size_t tile_size = std::min(IND_TILE, n_ind_ - ind_tile);

      // decompress the standardised genotype matrix in manageable tiles
      geno_buf.decompress(
        ind_tile,
        ind_tile + tile_size,
        pheno_loc,
        gen_tile_,
        true);

      // update the genetic component buffer in-place
      pheno_gen_buf.col(pheno).segment(ind_tile, tile_size).noalias() =
          std::sqrt(h2_gen_[pheno]) *
          gen_tile_.topRows(tile_size).leftCols(n_causal) * pheno_effects;
    }
  }
}

void ScorePhenotypes::scoreEnvironmental(State& state) {
  auto pheno_env_buf = state.pheno(ComponentType::ENVIRONMENTAL);
  rng::NormalPolar::fill(
      state.pheno(ComponentType::ENVIRONMENTAL).data(), n_ind_ * n_pheno_);
  pheno_env_buf = pheno_env_buf * env_chol_.transpose();

  for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
    pheno_env_buf.col(pheno).noalias() =
        std::sqrt(h2_env_[pheno]) * pheno_env_buf.col(pheno);
}

void ScorePhenotypes::scoreTotal(State& state) {
  state.pheno(ComponentType::TOTAL).noalias() =
      state.pheno(ComponentType::GENETIC) +
      state.pheno(ComponentType::ENVIRONMENTAL);
}

void ScorePhenotypes::operator()(State& state) {
  scoreGenetic(state);
  scoreEnvironmental(state);
  scoreTotal(state);
}

}  // namespace amsim::phenome
