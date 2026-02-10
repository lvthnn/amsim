#include <amsim/initialise.h>
#include <amsim/output/estimator.h>
#include <amsim/params.h>
#include <amsim/preprocess.h>
#include <amsim/state.h>
#include <amsim/transform.h>

#include <filesystem>
#include <iostream>

namespace amsim {

void setup_output(const std::filesystem::path& out_dir) {
  if (std::filesystem::exists(out_dir))
    throw std::runtime_error(
        "directory " + out_dir.string() + " already exists!");
  std::filesystem::create_directory(out_dir);
}

void simulation_preprocess(Params& params) {
  preprocess::OptimisePhenotypeArchitecture opt(params);

  opt();
}

// runs a single-threaded simulation
// this should be used by multithreaded controller
void simulation_run(
    State& state,
    Params& params,
    const Estimators& estimators,
    std::size_t n_gen,
    std::optional<std::size_t> rep_id) {
  // set up output directory
  setup_output(params.sim.out_dir);

  rng::set_seed(rng::auto_seed(params.sim.rng_seed));

  simulation_preprocess(params);

  // founder haplotype initialiser
  genome::HaplotypeGeneratorIID haplo(params);

  // set up required transformers
  phenome::ScorePhenotypes score(params);
  mating::AssortativeMating mate(params);
  genome::UpdateGenome update(params);

  // estimator super-transformer
  ComputeEstimates estimate(params, estimators, rep_id);

  // generate the initial state
  std::cout << "generating haplotypes\n";
  haplo.generate_haplotypes(state.geno);

  // stream headers to estimator files
  estimate.headers();

  while (state.gen < n_gen) {
    state.geno.compute_mafs();
    state.geno.compute_stats();

    // score phenotypes using operator()
    score(state);
    state.pheno.compute_stats();

    // match mates
    mate(state);

    // TODO(karihlynsson): Produce subviews / sample the population
    //     for (each view of the population)
    //       std::size_t view_mask = view->sample(state);
    //       estimate(state, view);

    // compute metrics for everything
    estimate(state);

    // Update the genome
    state.geno.transpose();
    update(state);
    state.geno.transpose();

    ++state.gen;
  }
}

// void simulations_run(
//   State& state,
//
// )

}  // namespace amsim
