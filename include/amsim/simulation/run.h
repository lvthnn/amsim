#pragma once

#include <amsim/core.h>
#include <amsim/estimate.h>
#include <amsim/init.h>
#include <amsim/sample.h>
#include <amsim/transform.h>

#include <atomic>
#include <iostream>
#include <stdexcept>
#include <thread>

namespace amsim {

namespace details {

inline std::uint64_t shuffle_seed(std::uint64_t rng_seed, std::size_t rep_id) {
  constexpr std::uint64_t Phi = 0x9E3779B97F4A7C15ULL;
  return rng_seed + (Phi * rep_id);
}

inline Params preprocess_simulation(const Simulation& simulation) {
  Params params = build_params(simulation);
  OptimisePhenotypeArchitecture opt(params);
  opt();
  std::cout << opt.expected() << "\n";

  return params;
}

inline void setup_outdir(const Simulation& simulation) {
  if (std::filesystem::exists(simulation.output_dir))
    throw std::runtime_error(
        "Directory " + simulation.output_dir.string() + "already exists!");
  std::filesystem::create_directory(simulation.output_dir);
}

inline void setup_log(const Simulation& simulation) {
  if (simulation.log_file)
    Log::file(simulation.output_dir / "amsim.log", simulation.log_level);
  else
    Log::stream(std::cout, simulation.log_level);
}

inline std::filesystem::path setup_replicate(
    const Simulation& simulation, std::size_t rep_id) {
  return simulation.output_dir / std::format("rep_{:03}", rep_id + 1);
}

inline void setup_writer(
    const Simulation& simulation, std::size_t n_replicates) {
  Writer::get_instance(
      simulation.output_dir, simulation.n_generations, n_replicates);
}
}  // namespace details

inline void run_simulation(
    const Params& params,
    const std::vector<PopulationEstimator>& estimators,
    const std::vector<SampleSpec>& samples) {
  rng::set_seed(params.sim.rng_seed);

  if (!std::filesystem::exists(params.sim.out_dir))
    std::filesystem::create_directory(params.sim.out_dir);

  State state = build_state(params);

  ComputePopulationEstimates population_est(params, estimators);
  ComputeSampleEstimates sample_est(params, samples);

  HaplotypeGeneratorIID haplo(params);
  RandomMating founder_mate(params);
  ScorePhenotypes score(params);
  AssortativeMating mate(params);
  UpdateGenome update(params);

  // initial state generation
  haplo(state);
  state.geno().compute_mafs();
  state.geno().compute_stats();

  score(state);
  state.pheno().compute_stats();
  founder_mate(state);
  state.transpose();
  update(state);
  state.transpose();
  state.advance();

  // run the core simulation loop
  while (state.gen <= params.sim.n_gens) {
    Log::debug("Simulating generation " + std::to_string(state.gen));
    state.geno().compute_mafs();
    state.geno().compute_stats();

    // score phenotypes using operator()
    score(state);
    state.pheno().compute_stats();

    // match mates
    mate(state);

    // sample and estimate subpopulations
    sample_est(state);

    // compute population estimates
    population_est(state);

    // update the genome — individual-major
    state.transpose();
    update(state);
    state.transpose();

    state.advance();
  }
}

inline void run_simulations(
    const Simulation& simulation,
    std::size_t n_replicates,
    std::size_t n_threads) {
  // set up thread pool for parallel simulation
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> pool;
  pool.reserve(n_threads);

  // set up random seed
  std::uint64_t seed = rng::auto_seed(simulation.random_seed);
  rng::set_seed(seed);

  // set up the base and replicate directories
  details::setup_outdir(simulation);

  // preprocess parameters
  Params params = details::preprocess_simulation(simulation);

  details::setup_log(simulation);
  details::setup_writer(simulation, n_replicates);

  for (std::size_t thread = 0; thread < n_threads; ++thread) {
    pool.emplace_back([&, thread]() {
      while (true) {
        Params thread_params = params;

        Log::debug("setting up thread " + std::to_string(thread));

        auto rep_id = next.fetch_add(1);
        if (rep_id >= n_replicates) return;

        thread_params.sim.out_dir =
            details::setup_replicate(simulation, rep_id);
        thread_params.sim.rng_seed = details::shuffle_seed(seed, rep_id);
        thread_params.sim.rep_id = rep_id;

        run_simulation(
            thread_params, simulation.estimators, simulation.samples);
      }
    });
  }

  for (auto& t : pool) t.join();
}

}  // namespace amsim
