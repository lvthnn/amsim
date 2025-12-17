#include <amsim/log_level.h>
#include <amsim/logger.h>
#include <amsim/logger_timer.h>
#include <amsim/metric.h>
#include <amsim/simulation.h>
#include <amsim/simulation_config.h>
#include <sys/resource.h>

#include <atomic>
#include <filesystem>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>

namespace amsim {

void resolve_rlimit(const std::size_t required) {
  auto rl_req = static_cast<rlim_t>(required);
  struct rlimit rl;

  if (getrlimit(RLIMIT_NOFILE, &rl) != 0) perror("getrlimit");

  if (rl_req > rl.rlim_max)
    throw std::invalid_argument(
        "required resource limit exceeds hard limit; consider lowering "
        "number of threads or the number of metrics");

  if (rl_req > rl.rlim_cur) {
    rl.rlim_cur = rl_req;
    if (setrlimit(RLIMIT_NOFILE, &rl) != 0) perror("setrlimit");
  }
}

Simulation::Simulation(
    const SimulationConfig& config,
    std::optional<std::filesystem::path> out_dir_,
    std::optional<std::uint64_t> rng_seed_)
    : n_gen(config.n_gen),
      n_ind(config.n_ind),
      n_loc(config.n_loc),
      n_pheno(config.n_pheno),
      pheno_names(config.v_name),
      out_dir(out_dir_ ? *out_dir_ : config.out_dir),
      genome_(n_ind, n_loc, config.v_mut, config.v_rec, config.v_maf),
      arch_(
          n_pheno,
          n_loc,
          config.v_n_loc,
          config.v_name,
          config.v_h2_gen,
          config.v_h2_env,
          config.v_h2_vert,
          config.v_rvert_pat,
          config.v_rvert_env,
          config.gen_cor,
          config.env_cor),
      buf_(n_ind, n_pheno, config.require_lat),
      phenotypes_([&]() {
        PhenotypeList phenotypes;
        phenotypes.reserve(n_pheno);
        for (std::size_t pheno = 0; pheno < n_pheno; ++pheno) {
          phenotypes.emplace_back(buf_, arch_, pheno_names[pheno], pheno);
        }
        return phenotypes;
      }()),
      model_(
          phenotypes_,
          config.mate_cor,
          n_ind / 2,
          config.n_itr,
          config.temp_init,
          config.temp_decay,
          config.tol_inf),
      ctx_(genome_, arch_, buf_, phenotypes_, model_),
      metrics_([&]() {
        std::vector<Metric> metrics;
        metrics.reserve(config.specs.size());
        for (const MetricSpec& spec : config.specs)
          metrics.push_back(spec.setup(ctx_));
        return metrics;
      }()) {
  if (!std::filesystem::exists(out_dir))
    std::filesystem::create_directory(out_dir);
  rng::set_seed(rng::auto_seed(rng_seed_ ? *rng_seed_ : config.rng_seed));
}

void Simulation::stream(std::size_t gen) {
  if (gen == 0) {
    streams_.reserve(metrics_.size());

    for (auto& metric : metrics_) {
      auto stream =
          std::make_unique<std::ofstream>(out_dir / (metric.name + ".tsv"));
      if (!stream->is_open()) {
        throw std::runtime_error(
            "could not open metric stream: " + metric.name);
      }
      *stream << metric.header() << "\n";
      streams_.emplace_back(std::move(stream));
    }
  }

  for (std::size_t metric = 0; metric < metrics_.size(); ++metric) {
    if (!streams_[metric] || !streams_[metric]->is_open()) {
      streams_[metric] = std::make_unique<std::ofstream>(
          out_dir / (metrics_[metric].name + ".tsv"));
      if (!streams_[metric]->is_open())
        throw std::runtime_error(
            "could not open metric file: " + metrics_[metric].name);
      *streams_[metric] << metrics_[metric].header() << "\n";
    }
    *streams_[metric] << (gen + 1) << "\t" << metrics_[metric].stream(ctx_)
                      << "\n";
  }
}

void Simulation::run() {
  LOG_DEBUG("Setting up simulation");
  std::vector<std::size_t> sib_matching(ctx_.n_ind / 2);
  std::iota(sib_matching.begin(), sib_matching.end(), ctx_.n_ind / 2);

  LOG_DEBUG("Generating haplotypes");
  genome_.generate_haplotypes();

  for (std::size_t gen = 0; gen < n_gen; ++gen) {
    LOG_DEBUG("Simulating generation " + std::to_string(gen));
    genome_.compute_mafs();
    genome_.compute_stats();

    arch_.gen_env(buf_(ComponentType::ENVIRONMENTAL), ctx_.n_ind);

    for (Phenotype& pheno : phenotypes_) {
      pheno.score(genome_);
      if (gen == 0) {
        amsim::PhenoArch::gen_vert(
            pheno(ComponentType::VERTICAL), n_ind, pheno.h2_vert());
      }
      pheno.score_tot();
      pheno.compute_stats();
    }

    if (buf_.has_lat()) buf_.score_latent(model_.cor_U, model_.cor_VT);

    model_.init_state();
    model_.update(phenotypes_);
    std::vector<std::size_t> opt_matching = model_.match();

    stream(gen);

    genome_.transpose();
    genome_.update(opt_matching, arch_, buf_);
    genome_.transpose();
  }
}

std::uint64_t shuffle_seed(std::uint64_t rng_seed, std::size_t rep_id) {
  constexpr std::uint64_t PHI = 0x9E3779B97F4A7C15ULL;
  return rng_seed + (PHI * rep_id);
}

void run_simulation(
    const SimulationConfig& config,
    std::optional<std::filesystem::path> out_dir_,
    bool log_file,
    LogLevel log_level) {
  Simulation sim(config, out_dir_);

  if (log_file)
    LOG_FILE((out_dir_ ? *out_dir_ : config.out_dir), log_level);
  else
    LOG_STREAM(std::cout, log_level);

  sim.run();
}

void run_simulations(
    const SimulationConfig& config,
    std::size_t n_replicates,
    std::size_t n_threads,
    bool summarise,
    bool log_file,
    LogLevel log_level) {
  // ensure the base directory exists
  if (!std::filesystem::exists(config.out_dir))
    std::filesystem::create_directory(config.out_dir);

  // set up the logger instance
  if (log_file)
    LOG_FILE(config.out_dir / "amsim.log", log_level);
  else
    LOG_STREAM(std::cout, log_level);

  LOG_INFO("Starting simulation");
  LOG_DEBUG("Using random seed " + std::to_string(config.rng_seed));

  // run multithreaded replicate simulations
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> pool;
  pool.reserve(n_threads);

  std::size_t n_metrics = config.specs.size();

  // set rlimit or warn user if process hard limit exceeded
  resolve_rlimit((n_metrics * n_threads) + 100);

  for (std::size_t ii = 0; ii < n_threads; ++ii) {
    pool.emplace_back([&]() {
      while (true) {
        std::size_t rep_id = next.fetch_add(1);
        if (rep_id >= n_replicates) return;

        // adjust the output directory path for the replicate
        std::filesystem::path rep_dir =
            config.out_dir / std::format("rep_{:03}", rep_id);

        // adjust replicate id rng seed
        std::uint64_t rep_seed = shuffle_seed(config.rng_seed, rep_id);

        Simulation rep(config, rep_dir, rep_seed);

        rep.run();
      }
    });
  }

  for (std::thread& thread : pool) thread.join();

  if (summarise) {
    LOG_INFO("Summarising results");
    SimulationResults results(config.out_dir);
    results.summarise();
    results.save();
  }

  LOG_INFO("Done!");
}

}  // namespace amsim
