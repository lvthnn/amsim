// This example shows the setup and running of a toy simulation configuration
#include <amsim/simulation.h>
#include <amsim/core.h>
#include <Eigen/Dense>

int main() {
  amsim::utils::check_plink2();
  amsim::utils::check_gcta64();

  Eigen::MatrixXd genetic_component_cor(2, 2);
  genetic_component_cor << 1.0, 0.5, 0.5, 1.0;

  Eigen::MatrixXd environmental_component_cor(2, 2);
  environmental_component_cor << 1.0, 0.5, 0.5, 1.0;

  Eigen::MatrixXd mate_cor(2, 2);
  mate_cor << 0.3, 0.5, 0.2, 0.4;

  auto simulation = amsim::Simulation{
      .n_individuals = 10000,
      .genome =
          amsim::Genome{
              .n_loci = 5000, .v_rec = 0.5, .v_maf = 0.5, .v_mut = 0.0},
      .phenotypes =
          {amsim::Phenotype{
               .name = "height",
               .n_causal_loci = 2000,
               .h2_genetic = 0.5,
               .h2_environmental = 0.25,
               .h2_nurture = 0.25},
           amsim::Phenotype{
               .name = "weight",
               .n_causal_loci = 2000,
               .h2_genetic = 0.75,
               .h2_environmental = 0.125,
               .h2_nurture = 0.125}},
      .genetic_component_cor = genetic_component_cor,
      .environmental_component_cor = environmental_component_cor,
      .mating = {.mate_cor = mate_cor},
      .estimators =
          {amsim::PopulationHeritability(),
           amsim::PopulationMateCor(),
           amsim::PopulationComponentVar(amsim::Component::Genetic),
           amsim::PopulationComponentVar(amsim::Component::Nurture),
           amsim::PopulationComponentVar(amsim::Component::Total),
           amsim::PopulationComponentCor(amsim::Component::Genetic),
           amsim::PopulationComponentCor(amsim::Component::Environmental),
           amsim::PopulationComponentCor(
               amsim::Component::Genetic, amsim::Component::Nurture)},
      .samples =
          {amsim::Sample<amsim::Proband::Individual>{
               .name = "the_broken_faucet",
               .n_probands = 200,
               .weighting = amsim::Uniform(),
               .estimators =
                   {amsim::SampleMeanEstimator<amsim::Proband::Individual>(),
                    amsim::SampleVarEstimator<amsim::Proband::Individual>()}},
           amsim::Sample<amsim::Proband::Family>{
               .name = "jonas_hamburger_study",
               .n_probands = 150,
               .weighting = amsim::Uniform(),
               .estimators =
                   {amsim::SampleMeanEstimator<amsim::Proband::Family>(),
                    amsim::SampleVarEstimator<amsim::Proband::Family>(),
                    amsim::SampleMateCorEstimator<amsim::Proband::Family>()}},
           amsim::Sample<amsim::Proband::Mate>{
               .name = "the_beautiful_people",
               .n_probands = 150,
               .of = amsim::Mate::All,
               .weighting = amsim::Uniform(),
               .estimators =
                   {amsim::SampleMateCorEstimator<amsim::Proband::Mate>()}}},
      .n_generations = 1,
      .output_dir = "amsim_demo"};

  amsim::run_simulations(simulation, 1, 1);
}
