#include <amsim/core.h>
#include <amsim/init.h>
#include <amsim/estimate.h>
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
  OptimisePhenotypeArchitecture opt(params);
  opt();
}

// runs a single-threaded simulation
// this should be used by multithreaded controller
void simulation_run(
    const Simulation& simulation,
    std::size_t n_gen,
    std::optional<std::size_t> rep_id) {
  if (simulation.log_file)
    LOG_FILE(simulation.output_dir, simulation.log_level);
  else
    LOG_STREAM(std::cout, simulation.log_level);

  LOG_DEBUG("Building parameters and preprocessing...");

  // build parameters, preprocess, and build state
  Params params = build_params(simulation);
  rng::set_seed(rng::auto_seed(params.sim.rng_seed));
  simulation_preprocess(params);

  State state = build_state(params);

  // set up output directory
  setup_output(params.sim.out_dir);

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
  HaplotypeGeneratorIID haplo(params);
  RandomMating random_mate(params);
  ScorePhenotypes score(params);
  AssortativeMating mate(params);
  UpdateGenome update(params);

  LOG_DEBUG("Generating initial state...");

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

  LOG_DEBUG("Finished generating initial state...");

  // run the core simulation loop
  while (state.gen <= n_gen) {
    LOG_DEBUG("Simulating generation " + std::to_string(state.gen));
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
