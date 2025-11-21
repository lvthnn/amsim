#include <amsim/component_type.h>
#include <amsim/mating.h>
#include <amsim/metricspec.h>
#include <amsim/simulation.h>
#include <amsim/simulation_config.h>
#include <amsim/logger.h>

#include <cstddef>
#include <vector>

int main() {
  std::vector<double> v_maf(500, 0.5);
  std::vector<double> v_rec(500, 0.5);
  std::vector<double> v_mut(500, 0.0);

  amsim::SimulationConfig config;

  config.simulation(10, 256000, "amsim_multithread", 1234568ull);

  config.genome(500, v_maf, v_rec, v_mut);

  config.phenome(
      2,
      {"height", "weight"},
      {250, 250},
      {0.5, 0.5},
      {0.5, 0.5},
      {0.0, 0.0},
      {1.0, 0.5, 0.5, 1.0},
      {1.0, 0.5, 0.5, 1.0});

  config.assortative_mating(std::vector<double>{0.4, 0.3, 0.2, 0.5}, 1e-7, 2e6);

  config.metrics(
      {amsim::pheno_h2(),
       amsim::pheno_comp_cor(amsim::ComponentType::GENETIC),
       amsim::pheno_comp_cor(amsim::ComponentType::ENVIRONMENTAL),
       amsim::pheno_comp_xcor(amsim::ComponentType::GENETIC),
       amsim::pheno_comp_xcor(amsim::ComponentType::TOTAL),
       amsim::pheno_latent_h2(),
       amsim::pheno_latent_comp_cor(amsim::ComponentType::GENETIC),
       amsim::pheno_latent_comp_xcor(amsim::ComponentType::GENETIC)});

  std::size_t n_replicates = 100;
  std::size_t n_threads = 10;

  amsim::run_simulations(
      config, n_replicates, n_threads, true, false, amsim::LogLevel::DEBUG);
}
