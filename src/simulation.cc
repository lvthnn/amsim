#include <amsim/initialise.h>
#include <amsim/output/estimator.h>
#include <amsim/params.h>
#include <amsim/preprocess.h>
#include <amsim/setup.h>
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
    const Simulation& simulation,
    std::size_t n_gen,
    std::optional<std::size_t> rep_id) {
  // build parameters, preprocess, and build state
  Params params = build_params(simulation);
  simulation_preprocess(params);
  State state = build_state(params);

  for (std::size_t pheno = 0; pheno < params.pheno.n_pheno; ++pheno) {
    std::cout << params.pheno.names[pheno] << "\n";
  }

  for (std::size_t pheno = 0; pheno < params.pheno.n_pheno; ++pheno) {
    std::cout << params.pheno.pheno_ids.at(params.pheno.names[pheno]) << "\n";
  }

  // set up output directory
  setup_output(params.sim.out_dir);
  rng::set_seed(rng::auto_seed(params.sim.rng_seed));

  ComputePopulationEstimates estimate(params, simulation.estimators, rep_id);

  // convert the specification objects into samplers
  std::vector<Sampler> samplers;
  samplers.reserve(simulation.samples.size());
  for (const auto& sample : simulation.samples) {
    samplers.push_back(
        std::visit(
            [&params](auto&& spec) { return Sampler(spec, params); }, sample));
  }


  // founder haplotype initialiser
  genome::HaplotypeGeneratorIID haplo(params);
  mating::RandomMating random_mate(params);
  phenome::ScorePhenotypes score(params);
  mating::AssortativeMating mate(params);
  genome::UpdateGenome update(params);

  // generate the initial state
  haplo.generate_haplotypes(state.geno());
  state.geno().compute_mafs();
  state.geno().compute_stats();

  score(state);
  state.pheno().compute_stats();
  random_mate(state);
  state.transpose();
  update(state);
  state.transpose();
  state.advance();

  // run the core simulation loop
  while (state.gen <= n_gen) {
    state.geno().compute_mafs();
    state.geno().compute_stats();

    // score phenotypes using operator()
    score(state);
    state.pheno().compute_stats();

    // match mates
    mate(state);

    // sample and estimate subpopulations
    for (auto& sampler : samplers)
      sampler(state);

    // compute population estimates
    estimate(state);

    // update the genome
    state.transpose();
    update(state);
    state.transpose();

    state.advance();
  }
}

}  // namespace amsim
