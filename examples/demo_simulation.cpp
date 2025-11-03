#include <amsim/component_type.h>
#include <amsim/mating.h>
#include <amsim/metricspec.h>
#include <amsim/simulation_config.h>
#include <amsim/simulation.h>

#include <cstddef>
#include <vector>
#include <iostream>

int main() {
  std::cout << "Alloc vectors..." << std::endl;
  std::vector<double> v_maf(1000, 0.5);
  std::vector<double> v_rec(1000, 0.5);
  std::vector<double> v_mut(1000, 1e-8);

  std::cout << "Create config..." << std::endl;
  amsim::SimulationConfig config;

  std::cout << "Config: simulation()" << std::endl;
  config.simulation(10, 5000, "amsim_multithread", 1234568ull);

  std::cout << "Config: genome()" << std::endl;
  config.genome(1000, v_maf, v_rec, v_mut);

  std::cout << "Config: phenome()" << std::endl;
  config.phenome(
      2,
      {"height", "weight"},
      {500, 500},
      {0.5, 0.5},
      {0.5, 0.5},
      {0.0, 0.0},
      {1.0, 0.5, 0.5, 1.0},
      {1.0, 0.5, 0.5, 1.0});

  std::cout << "Config: mating()" << std::endl;
  config.mating(
      amsim::MatingType::ASSORTATIVE,
      2e6,
      0.5,
      0.9999,
      std::vector<double>{0.4, 0.3, 0.2, 0.5});

  std::cout << "Config: metrics()" << std::endl;
  config.metrics(
      {amsim::pheno_h2(),
       amsim::pheno_comp_cor(amsim::ComponentType::GENETIC),
       amsim::pheno_comp_cor(amsim::ComponentType::ENVIRONMENTAL),
       amsim::pheno_comp_xcor(amsim::ComponentType::GENETIC),
       amsim::pheno_comp_xcor(amsim::ComponentType::TOTAL),
       amsim::pheno_latent_h2(),
       amsim::pheno_latent_comp_cor(amsim::ComponentType::GENETIC),
       amsim::pheno_latent_comp_xcor(amsim::ComponentType::GENETIC)});

  std::cout << "Run simulations..." << std::endl;
  std::size_t n_replicates = 10;
  std::size_t n_threads = 5;

  amsim::run_simulations(config, n_replicates, n_threads);

  std::cout << "Done." << std::endl;

}
