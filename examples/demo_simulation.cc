// This example shows the setup and running of a toy simulation configuration
#include <amsim/output/estimator.h>
#include <amsim/simulation.h>
#include <amsim/setup.h>
#include <amsim/state.h>

#include <Eigen/Dense>

int main() {
  Eigen::MatrixXd genetic_component_cor(2, 2);
  genetic_component_cor << 1.0, 0.5, 0.5, 1.0;

  Eigen::MatrixXd environmental_component_cor(2, 2);
  environmental_component_cor << 1.0, 0.5, 0.5, 1.0;

  Eigen::MatrixXd mate_cor(2, 2);
  mate_cor << 0.2, 0.4, 0.3, 0.5;

  auto simulation = amsim::Simulation{
    .n_individuals = 4000,
    .genome = amsim::Genome{
      .n_loci = 4000,
      .v_rec = 0.5,
      .v_maf = 0.5,
      .v_mut = 0.0
    },
    .phenotypes = {
      amsim::Phenotype{.name = "height", .n_causal_loci = 2000},
      amsim::Phenotype{.name = "weight", .n_causal_loci = 2000}
    },
    .genetic_component_cor = genetic_component_cor,
    .environmental_component_cor = environmental_component_cor,
    .mating = {.mate_cor = mate_cor},
    .estimators = {
      amsim::Heritability(),
      amsim::MateCor(),
      amsim::ComponentCor(amsim::phenome::Component::Genetic),
      amsim::ComponentCor(amsim::phenome::Component::Environmental)
    },
    .output_dir = "amsim_demo"
  };

  // in the future, we might consider adding
  //   auto pipeline = amsim::build_pipeline(simulation)
  auto params = amsim::build_params(simulation);
  auto state = amsim::build_state(params);

  amsim::simulation_run(state, params, simulation.estimators, 25);

  // run multiple simulations
  // amsim::run_simulations(simulation, 100, 10);
}
