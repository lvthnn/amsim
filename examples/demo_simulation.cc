// This example shows the setup and running of a toy simulation configuration
#include <amsim/data/genome.h>
#include <amsim/data/mating.h>
#include <amsim/data/phenome.h>
#include <amsim/output/estimator.h>
#include <amsim/params.h>
#include <amsim/simulation.h>
#include <amsim/state.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cstddef>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <iostream>

int main() {
  // Simulation settings
  constexpr std::size_t N_IND = 64000;
  constexpr std::size_t N_LOC = 2000;
  constexpr std::size_t N_GEN = 10;
  constexpr std::size_t N_PHENO = 2;
  std::string output_dir = "amsim_multithread";

  // Genome parameters: MAF, recombination, and mutation rates
  Eigen::VectorXd v_maf = Eigen::VectorXd::Constant(N_LOC, 0.5);
  Eigen::VectorXd v_rec = Eigen::VectorXd::Constant(N_LOC, 0.5);
  Eigen::VectorXd v_mut = Eigen::VectorXd::Constant(N_LOC, 0.0);

  amsim::genome::GenomeParams geno_params(N_IND, N_LOC, v_maf, v_rec, v_mut);

  // Phenome parameters: two phenotypes ("height" and "weight")
  std::vector<std::string> pheno_names = {"height", "weight"};
  std::vector<std::size_t> n_locs_per_pheno = {1000, 1000};

  Eigen::VectorXd h2_gen(N_PHENO);
  h2_gen << 0.5, 0.5;

  Eigen::VectorXd h2_env(N_PHENO);
  h2_env << 0.5, 0.5;

  Eigen::VectorXd h2_nur(N_PHENO);
  h2_nur << 0.0, 0.0;

  // Genetic and environmental correlation matrices
  Eigen::MatrixXd gen_cor(N_PHENO, N_PHENO);
  gen_cor << 1.0, 0.5, 0.5, 1.0;

  Eigen::MatrixXd env_cor(N_PHENO, N_PHENO);
  env_cor << 1.0, 0.5, 0.5, 1.0;

  amsim::phenome::PhenomeParams pheno_params(
      N_PHENO,
      pheno_names,
      n_locs_per_pheno,
      h2_gen,
      h2_env,
      h2_nur,
      gen_cor,
      env_cor);

  // Assign random causal loci and uniform effect sizes for each phenotype
  std::mt19937 rng(42);
  std::vector<std::size_t> all_loci(N_LOC);
  std::iota(all_loci.begin(), all_loci.end(), 0);

  for (std::size_t p = 0; p < N_PHENO; ++p) {
    // Shuffle and take the first n_locs_per_pheno[p] loci
    std::shuffle(all_loci.begin(), all_loci.end(), rng);
    std::vector<std::size_t> loci(all_loci.begin(),
                                  all_loci.begin() + n_locs_per_pheno[p]);
    std::ranges::sort(loci);
    pheno_params.causal_loc.push_back(std::move(loci));

    // Uniform effect sizes (1 / sqrt(n_loci) for unit variance)
    double effect = 1.0 / std::sqrt(static_cast<double>(n_locs_per_pheno[p]));
    pheno_params.pheno_effects.emplace_back(
        Eigen::VectorXd::Constant(n_locs_per_pheno[p], effect));
  }

  // Mating parameters: assortative mating correlation matrix
  Eigen::MatrixXd mate_cor(N_PHENO, N_PHENO);
  mate_cor << 0.4, 0.3, 0.2, 0.5;

  amsim::mating::MatingParams mate_params(mate_cor);

  // Simulation parameters
  amsim::SimulationParams sim_params;
  sim_params.n_reps = 100;
  sim_params.n_threads = 10;
  sim_params.rng_seed = 1234568ULL;
  sim_params.out_dir = output_dir;
  sim_params.log_level = "DEBUG";

  // Bundle all parameters
  amsim::Params params(geno_params, pheno_params, mate_params, sim_params);

  // Set up estimators (metrics)
  amsim::Estimators estimators = {
    amsim::PhenoHeritability(),
    amsim::PhenoComponentMean(),
    amsim::PhenoComponentVar(),
    amsim::MateCorrelation(),
    amsim::PhenoComponentCor(
      amsim::phenome::ComponentType::GENETIC,
      amsim::phenome::ComponentType::ENVIRONMENTAL),
    amsim::PhenoComponentCor(amsim::phenome::ComponentType::GENETIC)
  };

  std::cout << "number of estimators: " << std::to_string(estimators.size()) << "\n";

  // Create data buffers for the simulation state
  amsim::genome::GenoBuf geno_buf(N_IND, N_LOC, v_mut, v_rec, v_maf);
  amsim::phenome::PhenoBuf pheno_buf(N_IND, N_PHENO);
  amsim::mating::Matching matching;

  // Create the simulation state
  amsim::State state{
      .geno = geno_buf, .pheno = pheno_buf, .matching = matching};

  // Run the simulation
  amsim::simulation_run(state, params, estimators, N_GEN);
}
