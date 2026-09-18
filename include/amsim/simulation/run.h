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
  Log::info("Preprocessing simulation");

  Params params = buildParams(spec);
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
  Log::debug("Created simulation workspace directory {}", tmp.string());
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

inline std::string describeGenomeParam(
    const std::variant<double, File<Eigen::MatrixXd>, Distribution>& v) {
  if (std::holds_alternative<double>(v))
    return std::format("constant({})", std::get<double>(v));
  if (std::holds_alternative<Distribution>(v))
    return std::format("distribution({})", std::get<Distribution>(v).name);
  return std::format(
      "file({})", std::get<File<Eigen::MatrixXd>>(v).path.string());
}

inline std::string describeMatrixParam(
    const std::optional<std::variant<File<Eigen::MatrixXd>, Eigen::MatrixXd>>&
        v) {
  if (!v.has_value()) return "default";
  if (std::holds_alternative<Eigen::MatrixXd>(v.value())) {
    const auto& m = std::get<Eigen::MatrixXd>(v.value());
    return std::format("inline matrix ({}x{})", m.rows(), m.cols());
  }
  return std::format(
      "file({})", std::get<File<Eigen::MatrixXd>>(v.value()).path.string());
}

inline void logSimulationSpec(const SimulationSpec& spec) {
  Log::debug("Starting simulation with config:");
  Log::debug("  n_individuals={}", spec.n_individuals);
  Log::debug("  n_generations={}", spec.n_generations);
  Log::debug("  n_replicates={}", spec.n_replicates);
  Log::debug("  n_threads={}", spec.n_threads);
  Log::debug(
      "  random_seed={}",
      spec.random_seed.has_value() ? std::to_string(spec.random_seed.value())
                                   : std::string("none (random)"));
  Log::debug("  share_init_state={}", spec.share_init_state);
  Log::debug(
      "  pedigree_warmup={}, pedigree_max_depth={}",
      spec.pedigree_warmup,
      spec.pedigree_max_depth);
  Log::debug("  output_dir={}", spec.output_dir.string());
  Log::debug("  output_name={}", spec.output_name.value_or("(default)"));
  Log::debug(
      "  log_level={}, log_to_file={}",
      logLevelToString(spec.log_level),
      spec.log_to_file);

  Log::debug("  genome.n_loci={}", spec.genome.n_loci);
  Log::debug(
      "  genome.locus_freq={}", describeGenomeParam(spec.genome.locus_freq));
  Log::debug(
      "  genome.locus_rec={}", describeGenomeParam(spec.genome.locus_rec));
  Log::debug(
      "  genome.locus_mut={}", describeGenomeParam(spec.genome.locus_mut));

  Log::debug("  mating.type={}", spec.mating.type);
  Log::debug("  mating.tolerance={}", spec.mating.tolerance);
  Log::debug("  mating.max_iterations={}", spec.mating.max_iterations);
  Log::debug(
      "  mating.initial_temperature={}", spec.mating.initial_temperature);
  Log::debug("  mating.temperature_decay={}", spec.mating.temperature_decay);
  Log::debug("  mating.mate_cor={}", describeMatrixParam(spec.mating.mate_cor));

  Log::debug(
      "  genetic_component_cor={}",
      describeMatrixParam(spec.genetic_component_cor));
  Log::debug(
      "  environmental_component_cor={}",
      describeMatrixParam(spec.environmental_component_cor));

  Log::debug("  phenotypes ({} total):", spec.phenotypes.size());
  for (const auto& pheno : spec.phenotypes) {
    Log::debug(
        "    {}: n_causal_loci={}, var_genetic={}, var_environmental={}, "
        "var_vertical={}",
        pheno.name,
        pheno.n_causal_loci,
        pheno.var_genetic,
        pheno.var_environmental,
        pheno.var_vertical);
  }

  Log::debug("  population estimators ({} total):", spec.estimators.size());
  for (const auto& est : spec.estimators) Log::debug("    {}", est.name);

  Log::debug("  sample specs ({} total):", spec.sample_spec.size());
  for (const auto& sample : spec.sample_spec) {
    Log::debug(
        "    {}: proband_type={}, n_probands={}, estimators=[{}]",
        sample.name,
        sample.proband_type,
        sample.n_probands,
        boost::algorithm::join(sample.estimators, ", "));
  }

  Log::debug(
      "  sample estimator specs ({} total):", spec.sample_estimator_spec.size());
  for (const auto& est : spec.sample_estimator_spec)
    Log::debug("    {} (type={})", est.name, est.type);
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
    Log::info("Simulating panmictic generation {}", it + 1);
    score(state);
    state.pheno().computeStats();
    founder_mate(state);
    state.updatePedigree();
    state.transpose();
    update(state);
    state.transpose();
    state.advance();
  }

  state.gen = 0;

  Log::info("Finished simulating founder generation(s)");

  if (params.global.post_init_seed.has_value()) {
    rng::setSeed(params.global.post_init_seed.value());
    Log::debug(
        "Setting post init seed from shared state to {}",
        params.global.post_init_seed.value());
  }

  // run the core simulation loop
  while (state.gen <= params.global.n_gens) {
    Log::info("Simulating generation {}", state.gen);
    state.geno().computeLocusFreqs();
    state.geno().computeLocusStats();

    // score phenotypes
    score(state);
    state.pheno().computeStats();

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

  details::logSimulationSpec(spec);

  // set up random seed
  std::uint64_t seed = rng::seedOrRandom(spec.random_seed);
  rng::setSeed(seed);

  // preprocess parameters
  Params params = details::preprocessSimulation(spec);
  params.global.rng_seed = seed;

  // validate output directory
  if (!std::filesystem::exists(spec.output_dir)) {
    try {
      std::filesystem::create_directory(spec.output_dir);
    } catch (std::exception& e) {
      Log::debug("Could not create output directory: {}", e.what());
      exit(EXIT_FAILURE);
    }
  }

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

  if (spec.delete_tmp)
    std::filesystem::remove_all(tmp_dir);
}

}  // namespace amsim
