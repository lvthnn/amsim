#include <cstddef>
#include <cstdint>
#include <string>

#include <amsim/simulation_builder.h>
#include <amsim/simulation.h>

namespace amsim {
  SimulationBuilder::SimulationBuilder()
    : status_() { }

  void SimulationBuilder::simulation(std::size_t n_gen, std::size_t n_ind,
                                     std::string out_path,
                                     std::uint64_t rng_seed) {
    n_gen_    = n_gen; 
    n_ind_    = n_ind;
    out_path_ = out_path;
    rng_seed_ = rng_seed;
  }

  void SimulationBuilder::genome(std::size_t n_loc, std::vector<double> v_maf,
                                 std::vector<double> v_rec,
                                 std::vector<double> v_mut) {
    n_loc_ = n_loc;
    v_maf_ = std::move(v_maf);
    v_rec_ = std::move(v_rec);
    v_mut_ = std::move(v_mut);
  }

  void SimulationBuilder::phenome(std::size_t n_pheno,
                                  std::vector<std::string> v_name,
                                  std::vector<double> v_h2_gen,
                                  std::vector<double> v_h2_env,
                                  std::vector<double> v_h2_vert,
                                  std::vector<double> gen_cor,
                                  std::vector<double> env_cor) {
    n_pheno_   = n_pheno;
    v_name_    = std::move(v_name);
    v_h2_gen_  = std::move(v_h2_gen);
    v_h2_env_  = std::move(v_h2_env);
    v_h2_vert_ = std::move(v_h2_vert);
    gen_cor_   = std::move(gen_cor);
    env_cor_   = std::move(env_cor);
  }

  void SimulationBuilder::mating(MatingType type,
                                 std::optional<std::size_t> n_itr,
                                 std::optional<double> tmp_init,
                                 std::optional<double> tmp_decay,
                                 std::optional<std::vector<double>> mate_cor) {
    type_ = type;
    if (type == MatingType::ASSORTATIVE) {

    }
  }

  Simulation SimulationBuilder::build() {

  }
}
