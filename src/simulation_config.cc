#include <amsim/logger.h>
#include <amsim/simulation_config.h>
#include <amsim/utils.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <stdexcept>

namespace amsim {

SimulationConfig& SimulationConfig::simulation(
    std::size_t n_gen_,
    std::size_t n_ind_,
    const std::string& out_dir_,
    std::optional<std::uint64_t> rng_seed_) {
  if (n_ind_ % 2 != 0) {
    n_ind_ += 1;
    LOG_INFO(
        "Rounding up population size to even number; " +
        std::to_string(n_ind_ - 1) + " -> " + std::to_string(n_ind));
  }
  n_gen = n_gen_;
  n_ind = n_ind_;
  out_dir = std::filesystem::path(out_dir_);

  if (rng_seed_)
    rng_seed = *rng_seed_;
  else
    rng_seed =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();

  return *this;
}

SimulationConfig& SimulationConfig::genome(
    std::size_t n_loc_,
    std::vector<double> v_maf_,
    std::vector<double> v_rec_,
    std::vector<double> v_mut_) {
  if (v_maf_.size() != n_loc_ || v_rec_.size() != n_loc_ ||
      v_mut_.size() != n_loc_)
    throw std::invalid_argument(
        "MAF, recombination, and mutation vectors must be of length equal to "
        "number of loci");

  utils::assert_probs(n_loc_, v_maf_.data(), 1);
  utils::assert_probs(n_loc_, v_rec_.data(), 1);
  utils::assert_probs(n_loc_, v_mut_.data(), 1);

  n_loc = n_loc_;
  v_maf = std::move(v_maf_);
  v_rec = std::move(v_rec_);
  v_mut = std::move(v_mut_);
  return *this;
}

SimulationConfig& SimulationConfig::phenome(
    std::size_t n_pheno_,
    const std::vector<std::string>& v_name_,
    std::vector<std::size_t> v_n_loc_,
    std::vector<double> v_h2_gen_,
    std::vector<double> v_h2_env_,
    std::vector<double> v_h2_vert_,
    std::vector<double> gen_cor_,
    std::vector<double> env_cor_) {
  if (v_name_.size() != n_pheno_)
    throw std::invalid_argument("must have equally many names as phenotypes");
  if (v_n_loc_.size() != n_pheno_)
    throw std::invalid_argument(
        "must have equally many number of loci per phenotype as phenotypes");
  if (v_h2_gen_.size() != n_pheno_)
    throw std::invalid_argument("must specify genetic h2 for all phenotypes");
  if (v_h2_env_.size() != n_pheno_)
    throw std::invalid_argument(
        "must specify environmental h2 for all phenotypes");
  if (v_h2_vert_.size() != n_pheno_)
    throw std::invalid_argument("must specify vertical h2 for all phenotypes");

  for (const std::size_t& nl : v_n_loc_)
    if (nl > n_loc)
      throw std::invalid_argument(
          "number of causal loci exceeds number of modelled loci");

  utils::assert_probs(n_pheno_, v_h2_gen_.data(), 1);
  utils::assert_probs(n_pheno_, v_h2_env_.data(), 1);
  utils::assert_probs(n_pheno_, v_h2_vert_.data(), 1);
  utils::assert_cor(n_pheno_, gen_cor_.data(), n_pheno_);
  utils::assert_cor(n_pheno_, env_cor_.data(), n_pheno_);

  for (std::size_t el = 0; el < n_pheno_; ++el)
    if (std::abs(v_h2_gen_[el] + v_h2_env_[el] + v_h2_vert_[el] - 1.0) > 1e-12)
      throw std::invalid_argument(
          "phenotype component variances must sum to one");

  n_pheno = n_pheno_;
  v_name = v_name_;
  v_n_loc = std::move(v_n_loc_);
  v_h2_gen = std::move(v_h2_gen_);
  v_h2_env = std::move(v_h2_env_);
  v_h2_vert = std::move(v_h2_vert_);
  gen_cor = std::move(gen_cor_);
  env_cor = std::move(env_cor_);
  return *this;
}

SimulationConfig& SimulationConfig::random_mating() {
  mating_type = MatingType::RANDOM;
  mate_cor = std::vector<double>(n_pheno * n_pheno, 0.0);
  tol_inf = 0.0;
  n_itr = 0;
  temp_init = 0.0;
  temp_decay = 0.0;
  return *this;
}

SimulationConfig& SimulationConfig::assortative_mating(
    std::vector<double> mate_cor_,
    std::optional<double> tol_inf_,
    std::optional<std::size_t> n_itr_,
    std::optional<double> temp_init_,
    std::optional<double> temp_decay_) {
  mating_type = MatingType::ASSORTATIVE;

  utils::assert_cross_cor(n_pheno, (mate_cor_).data(), n_pheno);
  mate_cor = std::move(mate_cor_);

  if (n_itr_) n_itr = *n_itr_;
  if (temp_init_) temp_init = *temp_init_;
  if (temp_decay_) temp_decay = *temp_decay_;
  if (tol_inf_) tol_inf = *tol_inf_;

  if (temp_init < 0.0)
    throw std::invalid_argument("temp_init must be positive");
  if (temp_decay < 0.0)
    throw std::invalid_argument("temp_decay must be positive");
  if (tol_inf < 0.0) throw std::invalid_argument("tol_inf must be positive");

  return *this;
}

SimulationConfig& SimulationConfig::metrics(std::vector<MetricSpec> specs_) {
  specs = std::move(specs_);
  require_lat = std::ranges::any_of(
      specs, [](const MetricSpec& s) { return s.require_lat; });
  return *this;
}

}  // namespace amsim
