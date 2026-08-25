// This file is part of amsim, copyright (C) 2025-2026 Kári Hlynsson.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <amsim/core.h>
#include <amsim/estimate.h>
#include <amsim/init.h>
#include <amsim/sample.h>
#include <amsim/transform.h>

#include <atomic>
#include <chrono>
#include <iostream>
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

inline std::filesystem::path setup_outdir(std::uint64_t seed) {
  auto now = std::chrono::system_clock::now().time_since_epoch();
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now).count();
  std::filesystem::path tmp =
      std::filesystem::temp_directory_path() /
      std::format("amsim_{}_{}", seed, ms);
  std::filesystem::create_directories(tmp);
  return tmp;
}

inline std::string results_filename(const Simulation& simulation) {
  return simulation.output_name.has_value()
      ? std::format("results_{}.h5", simulation.output_name.value())
      : "results.h5";
}

inline std::string log_filename(const Simulation& simulation) {
  return simulation.output_name.has_value()
      ? std::format("amsim_{}.log", simulation.output_name.value())
      : "amsim.log";
}

inline void setup_log(const Simulation& simulation) {
  if (simulation.log_to_file) {
    std::filesystem::path log_path =
        simulation.output_dir / log_filename(simulation);
    Log::file(log_path, simulation.log_level);
  } else {
    Log::stream(std::cout, simulation.log_level);
  }
}

inline std::filesystem::path setup_replicate(
    const std::filesystem::path& tmp_dir, std::size_t rep_id) {
  return tmp_dir / std::format("rep_{:03}", rep_id + 1);
}

inline void setup_writer(
    const Simulation& simulation, std::size_t n_replicates) {
  std::filesystem::path results_path =
      simulation.output_dir / results_filename(simulation);
  Writer::get_instance(results_path, simulation.n_generations, n_replicates);
}
}  // namespace details

inline void run_replicate(
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
  state.update_pedigree();
  state.transpose();
  update(state);
  state.transpose();
  state.advance();

  if (params.sim.post_init_seed.has_value())
    rng::set_seed(params.sim.post_init_seed.value());

  // run the core simulation loop
  while (state.gen <= params.sim.n_gens) {
    Log::info("Simulating generation " + std::to_string(state.gen));
    state.geno().compute_mafs();
    state.geno().compute_stats();

    // score phenotypes
    score(state);
    state.pheno().compute_stats();

    // match mates
    mate(state);
    state.update_pedigree();

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

inline void run_simulation(
    const Simulation& simulation,
    std::size_t n_replicates,
    std::size_t n_threads) {
  // set up thread pool for parallel simulation
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> pool;
  pool.reserve(n_threads);

  // set up the log
  details::setup_log(simulation);

  // set up random seed
  std::uint64_t seed = rng::auto_seed(simulation.random_seed);
  rng::set_seed(seed);

  // preprocess parameters
  Params params = details::preprocess_simulation(simulation);

  // validate output directory
  if (!std::filesystem::exists(simulation.output_dir))
    throw std::runtime_error(
        "output directory does not exist: " + simulation.output_dir.string());

  // set up scratch directory in tmp
  std::filesystem::path tmp_dir = details::setup_outdir(seed);

  details::setup_writer(simulation, n_replicates);

  Writer::get_instance().write_params(params);

  for (std::size_t thread = 0; thread < n_threads; ++thread) {
    pool.emplace_back([&, thread]() {
      while (true) {
        Params thread_params = params;

        auto rep_id = next.fetch_add(1);
        if (rep_id >= n_replicates) return;

        thread_params.sim.out_dir =
            details::setup_replicate(tmp_dir, rep_id);
        if (simulation.share_init_state) {
          thread_params.sim.rng_seed = seed;
          thread_params.sim.post_init_seed = details::shuffle_seed(seed, rep_id);
        } else {
          thread_params.sim.rng_seed = details::shuffle_seed(seed, rep_id);
        }
        thread_params.sim.rep_id = rep_id;

        run_replicate(
            thread_params, simulation.estimators, simulation.samples);
      }
    });
  }

  for (auto& t : pool) t.join();

  std::filesystem::remove_all(tmp_dir);
}

}  // namespace amsim
