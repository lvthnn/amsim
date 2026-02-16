// This example shows the setup and running of a toy simulation configuration
#include <amsim/init.h>
#include <amsim/core.h>
#include <amsim/estimate.h>
#include <amsim/simulation.h>

#include <Eigen/Dense>

int main() {
  Eigen::MatrixXd genetic_component_cor(2, 2);
  genetic_component_cor << 1.0, 0.5, 0.5, 1.0;

  Eigen::MatrixXd environmental_component_cor(2, 2);
  environmental_component_cor << 1.0, 0.5, 0.5, 1.0;

  Eigen::MatrixXd mate_cor(2, 2);
  mate_cor << 0.3, 0.5, 0.2, 0.4;

  auto simulation = amsim::Simulation{
    .n_individuals = 4000,
    .genome = amsim::Genome{
      .n_loci = 4000,
      .v_rec = 0.5,
      .v_maf = 0.5,
      .v_mut = 0.0
    },
    .phenotypes = {
      amsim::Phenotype{
        .name = "height",
        .n_causal_loci = 2000,
        .h2_genetic = 0.5,
        .h2_environmental = 0.25,
        .h2_nurture = 0.25
      },
      amsim::Phenotype{
        .name = "weight",
        .n_causal_loci = 2000,
        .h2_genetic = 0.75,
        .h2_environmental = 0.125,
        .h2_nurture = 0.125
      }
    },
    .genetic_component_cor = genetic_component_cor,
    .environmental_component_cor = environmental_component_cor,
    .mating = {.mate_cor = mate_cor},
    .estimators = {
      amsim::PopulationHeritability(),
      amsim::PopulationMateCor(),
      amsim::PopulationComponentVar(amsim::Component::Genetic),
      amsim::PopulationComponentVar(amsim::Component::Nurture),
      amsim::PopulationComponentVar(amsim::Component::Total),
      amsim::PopulationComponentCor(amsim::Component::Genetic),
      amsim::PopulationComponentCor(amsim::Component::Environmental),
      amsim::PopulationComponentCor(amsim::Component::Genetic,
                                    amsim::Component::Nurture)
    },
    .samples = {
      amsim::Sample<amsim::Proband::Individual>{
        .name = "random_sample_200",
        .n_probands = 200,
        .weighting = amsim::Uniform(),
        .estimators = {amsim::SampleMeanEstimator<amsim::Proband::Individual>()}
      },
      amsim::Sample<amsim::Proband::Individual>{
        .name = "participation_bias_200",
        .n_probands = 200,
        .on = std::vector<std::string>({"height", "weight"}),
        .weighting = amsim::Logistic((Eigen::VectorXd(2) << 1.0, -1.0).finished()),
        .estimators = {amsim::SampleMeanEstimator<amsim::Proband::Individual>()}
      }
    },
    .output_dir = "amsim_demo"
  };
  // in the future, we might consider adding
  //   auto pipeline = amsim::build_pipeline(simulation)
  amsim::simulation_run(simulation, 25);

  // run multiple simulations
  // amsim::run_simulations(simulation, 100, 10);
}
