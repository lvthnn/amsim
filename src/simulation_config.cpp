#include <amsim/simulation_config.h>

#include <filesystem>

namespace amsim {

SimulationConfig& SimulationConfig::simulation(
    std::size_t n_gen_,
    std::size_t n_ind_,
    std::string out_dir_,
    std::uint64_t rng_seed_) {
  n_gen = n_gen_;
  n_ind = n_ind_;
  out_dir = std::filesystem::path(out_dir_);
  rng_seed = rng_seed_;
  return *this;
}

SimulationConfig& SimulationConfig::genome(
    std::size_t n_loc_,
    std::vector<double> v_maf_,
    std::vector<double> v_rec_,
    std::vector<double> v_mut_) {
  n_loc = n_loc_;
  v_maf = std::move(v_maf_);
  v_rec = std::move(v_rec_);
  v_mut = std::move(v_mut_);
  return *this;
}

SimulationConfig& SimulationConfig::phenome(
    std::size_t n_pheno_,
    std::vector<std::string> v_name_,
    std::vector<std::size_t> v_n_loc_,
    std::vector<double> v_h2_gen_,
    std::vector<double> v_h2_env_,
    std::vector<double> v_h2_vert_,
    std::vector<double> gen_cor_,
    std::vector<double> env_cor_) {
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

SimulationConfig& SimulationConfig::mating(
    MatingType type_,
    std::optional<std::size_t> n_itr_,
    std::optional<double> temp_init_,
    std::optional<double> temp_decay_,
    std::optional<std::vector<double>> mate_cor_) {
  type = type_;

  if (type == MatingType::RANDOM) {
    n_itr = 0.0;
    temp_init = 0.0;
    temp_decay = 0.0;
    mate_cor = std::vector<double>(n_pheno * n_pheno, 0.0);
    return *this;
  }

  if (n_itr_) n_itr = *n_itr_;
  if (temp_init_) temp_init = *temp_init;
  if (temp_decay_) temp_decay = *temp_decay;
  if (!mate_cor) {
    throw std::runtime_error("need mate correlation for assortative model");
  }

  mate_cor_ = std::move(*mate_cor_);

  return *this;
}

SimulationConfig& SimulationConfig::metrics(std::vector<MetricSpec> specs_) {
  specs = std::move(specs_);
  require_lat =
      std::any_of(specs_.begin(), specs_.end(), [](const MetricSpec& s) {
        return s.require_lat;
      });
  return *this;
}

}  // namespace amsim
