// This example shows the setup and running of a toy simulation configuration
#include <amsim/component_type.h>
#include <amsim/mating.h>
#include <amsim/metricspec.h>
#include <amsim/simulation.h>
#include <amsim/simulation_config.h>
#include <amsim/logger.h>

#include <cstddef>
#include <string>
#include <vector>

int main() {
  std::string output_dir = "amsim_multithread";

  std::vector<double> v_maf(2000, 0.5);
  std::vector<double> v_rec(2000, 0.5);
  std::vector<double> v_mut(2000, 0.0);

  amsim::SimulationConfig config;

  config.simulation(10, 256000, output_dir, 1234568ULL);

  config.genome(2000, v_maf, v_rec, v_mut);

  config.phenome(
      2,
      {"height", "weight"},
      {1000, 1000},
      {0.5, 0.5},
      {0.5, 0.5},
      {0.0, 0.0},
      {1.0, 0.5, 0.5, 1.0},
      {1.0, 0.5, 0.5, 1.0});

  config.assortative_mating(std::vector<double>{0.4, 0.3, 0.2, 0.5}, 1e-7, 2e6);

  config.metrics(
      {amsim::pheno_h2(),
       amsim::pheno_comp_cor(amsim::ComponentType::GENETIC),
       amsim::loc_maf(),
       amsim::loc_var()});

  std::size_t n_replicates = 100;
  std::size_t n_threads = 10;

  amsim::run_simulations(
      config, n_replicates, n_threads, true, false, amsim::LogLevel::DEBUG);
}
