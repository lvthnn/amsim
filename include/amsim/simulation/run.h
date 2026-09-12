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

inline std::uint64_t shuffleSeed(std::uint64_t rng_seed, std::size_t rep_id) {
  constexpr std::uint64_t Phi = 0x9E3779B97F4A7C15ULL;
  return rng_seed + (Phi * rep_id);
}

inline Params preprocessSimulation(const SimulationSpec& spec) {
  Params params = build_params(spec);
  OptimisePhenotypeArchitecture opt(params);

  opt();

  return params;
}

inline std::filesystem::path outputSetup(std::uint64_t seed) {
  auto now = std::chrono::system_clock::now().time_since_epoch();
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now).count();
  std::filesystem::path tmp = std::filesystem::temp_directory_path() /
                              std::format("amsim_{}_{}", seed, ms);
  std::filesystem::create_directories(tmp);
  Log::debug(
      std::format("Created simulation workspace directory {}", tmp.string()));
  return tmp;
}

inline std::string resultsFilename(const SimulationSpec& spec) {
  return spec.output_name.has_value()
             ? std::format("results_{}.h5", spec.output_name.value())
             : "results.h5";
}

inline std::string logFilename(const SimulationSpec& spec) {
  return spec.output_name.has_value()
             ? std::format("amsim_{}.log", spec.output_name.value())
             : "amsim.log";
}

inline void logSetup(const SimulationSpec& spec) {
  if (spec.log_to_file) {
    std::filesystem::path log_path = spec.output_dir / logFilename(spec);
    Log::file(log_path, spec.log_level);
  } else {
    Log::stream(std::cout, spec.log_level);
  }
}

inline std::filesystem::path replicateSetup(
    const std::filesystem::path& tmp_dir, std::size_t rep_id) {
  return tmp_dir / std::format("rep_{:03}", rep_id + 1);
}

inline void writerSetup(const SimulationSpec& spec) {
  std::filesystem::path results_path = spec.output_dir / resultsFilename(spec);
  H5Writer::getInstance(results_path, spec.n_generations, spec.n_replicates);
}
}  // namespace details

inline void runReplicate(const Params& params) {
  rng::setSeed(params.global.rng_seed);

  if (!std::filesystem::exists(params.global.out_dir))
    std::filesystem::create_directory(params.global.out_dir);

  State state = buildState(params);

  ComputePopulationEstimates population_est(params);
  ComputeSampleEstimates sample_est(params);

  HaplotypeGeneratorIID haplo(params);
  RandomMating founder_mate(params);
  ScorePhenotypes score(params);
  AssortativeMating mate(params);
  UpdateGenome update(params);

  // initial state generation
  haplo(state);
  state.geno().computeLocusFreqs();
  state.geno().computeLocusStats();

  // if pedigree_warmup is true, then we simulate pedigree_max_depth
  // panmictic generations
  std::size_t n_it =
      params.global.pedigree_warmup ? params.global.pedigree_max_depth : 1;

  for (std::size_t it = 0; it < n_it; ++it) {
    Log::info(std::format("Simulating panmictic generation {}", it + 1));
    score(state);
    state.pheno().compute_stats();
    founder_mate(state);
    state.updatePedigree();
    state.transpose();
    update(state);
    state.transpose();
    state.advance();
  }

  Log::info("Finished simulating founder generation(s)");

  if (params.global.post_init_seed.has_value()) {
    rng::setSeed(params.global.post_init_seed.value());
    Log::debug(
        std::format(
            "Setting post init seed from shared state to {}",
            params.global.post_init_seed.value()));
  }

  // run the core simulation loop
  while (state.gen <= params.global.n_gens) {
    Log::info("Simulating generation " + std::to_string(state.gen));
    state.geno().computeLocusFreqs();
    state.geno().computeLocusStats();

    // score phenotypes
    score(state);
    state.pheno().compute_stats();

    // match mates
    mate(state);
    state.updatePedigree();

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

inline void runSimulation(const SimulationSpec& spec) {
  // set up thread pool for parallel simulation
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> pool;
  pool.reserve(spec.n_threads);

  // set up the log
  details::logSetup(spec);

  // set up random seed
  std::uint64_t seed = rng::seedOrRandom(spec.random_seed);
  rng::setSeed(seed);

  // preprocess parameters
  Params params = details::preprocessSimulation(spec);

  // validate output directory
  if (!std::filesystem::exists(spec.output_dir))
    throw std::runtime_error(
        "output directory does not exist: " + spec.output_dir.string());

  // set up scratch directory in tmp
  std::filesystem::path tmp_dir = details::outputSetup(seed);

  details::writerSetup(spec);

  H5Writer::getInstance().writeParams(params);

  for (std::size_t thread = 0; thread < spec.n_threads; ++thread) {
    pool.emplace_back([&, thread]() {
      while (true) {
        Params thread_params = params;

        auto rep_id = next.fetch_add(1);
        if (rep_id >= spec.n_replicates) return;

        thread_params.global.out_dir = details::replicateSetup(tmp_dir, rep_id);
        if (spec.share_init_state) {
          thread_params.global.rng_seed = seed;
          thread_params.global.post_init_seed =
              details::shuffleSeed(seed, rep_id);
        } else {
          thread_params.global.rng_seed = details::shuffleSeed(seed, rep_id);
        }
        thread_params.global.rep_id = rep_id;

        runReplicate(thread_params);
      }
    });
  }

  for (auto& t : pool) t.join();

  std::filesystem::remove_all(tmp_dir);
}

}  // namespace amsim
