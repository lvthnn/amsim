#include <cstddef>
#include <cstdint>
#include <string>
#include <iostream>
#include <filesystem>

#include <amsim/simulation_builder.h>
#include <amsim/simulation.h>

#include <amsim/rng.h>
#include <amsim/utils.h>
#include <amsim/genome.h>
#include <amsim/phenobuf.h>
#include <amsim/phenoarch.h>
#include <amsim/phenotype.h>
#include <amsim/mating.h>
#include <amsim/metric.h>

namespace amsim {
  SimulationBuilder::SimulationBuilder()
    : status_() { }

  SimulationBuilder& SimulationBuilder::simulation(std::size_t n_gen,
                                                   std::size_t n_ind,
                                                   std::string out_dir,
                                                   std::uint64_t rng_seed) {
    n_gen_    = n_gen; 
    n_ind_    = n_ind;
    out_dir_  = std::filesystem::path(out_dir);
    rng_seed_ = rng_seed;
    status_++;
    return *this;
  }

  SimulationBuilder& SimulationBuilder::genome(std::size_t n_loc,
                                               std::vector<double> v_maf,
                                               std::vector<double> v_rec,
                                               std::vector<double> v_mut) {
    if (status_ < SimulationStatus::GENOME)
      std::cerr << "Specify simulation parameters before configuring mating "
                   "model.\n";

    n_loc_ = n_loc;
    v_maf_ = std::move(v_maf);
    v_rec_ = std::move(v_rec);
    v_mut_ = std::move(v_mut);
    status_++;

    return *this;
  }

  SimulationBuilder& SimulationBuilder::phenome(std::size_t n_pheno,
                                                std::vector<std::string> v_name,
                                                std::vector<std::size_t> v_n_loc,
                                                std::vector<double> v_h2_gen,
                                                std::vector<double> v_h2_env,
                                                std::vector<double> v_h2_vert,
                                                std::vector<double> gen_cor,
                                                std::vector<double> env_cor) {
    if (status_ < SimulationStatus::PHENOME)
      std::cerr << "Specify simulation and genome parameters before configuring "
                   "mating model.\n";

    n_pheno_   = n_pheno;
    v_name_    = std::move(v_name);
    v_n_loc_   = std::move(v_n_loc);
    v_h2_gen_  = std::move(v_h2_gen);
    v_h2_env_  = std::move(v_h2_env);
    v_h2_vert_ = std::move(v_h2_vert);
    gen_cor_   = std::move(gen_cor);
    env_cor_   = std::move(env_cor);
    status_++;

    return *this;
  }

  SimulationBuilder& SimulationBuilder::mating(MatingType type,
                                               std::optional<std::size_t> n_itr,
                                               std::optional<double> tmp_init,
                                               std::optional<double> tmp_decay,
                                               std::optional<std::vector<double>> mate_cor) {
    if (status_ < SimulationStatus::MATING)
      std::cerr << "Specify simulation, genome, and phenome parameters before "
                   "configuring mating model.\n";

    if (type == MatingType::RANDOM) {
      n_itr_     = 0;
      tmp_init_  = 0; 
      tmp_decay_ = 0;
      mate_cor_  = std::vector<double>(n_pheno_ * n_pheno_, 0.0);
    }

    if (!n_itr || !tmp_init || !tmp_decay || !mate_cor)
      std::cerr << "Specify all of n_itr, tmp_init, tmp_decay, and mate_cor "
                   "for assortative mating models.\n";
    
    n_itr_     = *n_itr;
    tmp_init_  = *tmp_init;
    tmp_decay_ = *tmp_decay;
    mate_cor_  = *mate_cor;
    status_++;

    return *this;
  }

  SimulationBuilder& SimulationBuilder::metrics(std::vector<Metric> metrics) {
    metrics_ = std::move(metrics);
    return *this;
  }

  Simulation SimulationBuilder::build() {
    if (status_ < SimulationStatus::READY)
      std::cerr << "Specify simulation, genome, phenome, and mating parameters "
                   "before building simulation.\n";
		Simulation simulation(n_gen_, n_ind_, out_dir_, rng_seed_,
													n_loc_, v_maf_, v_rec_, v_mut_,
													n_pheno_, v_name_, v_n_loc_, v_h2_gen_, v_h2_env_, v_h2_vert_, gen_cor_, env_cor_,
													mate_cor_, n_itr_, tmp_init_, tmp_decay_,
													metrics_);
    return simulation;
  }
}
