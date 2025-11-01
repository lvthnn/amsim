#include <amsim/metric.h>
#include <amsim/simulation.h>
#include <amsim/simulation_config.h>

#include <atomic>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace amsim {

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
      rng_(
          rng::seed_xoshiro(
              rng::auto_seed(rng_seed_ ? *rng_seed_ : config.rng_seed))),
      genome_(n_ind, n_loc, config.v_mut, config.v_rec, config.v_maf, rng_),
      arch_(
          n_pheno,
          n_loc,
          config.v_n_loc,
          config.v_h2_gen,
          config.v_h2_env,
          config.gen_cor,
          config.env_cor,
          rng_),
      buf_(n_ind, n_pheno, config.require_lat),
      phenotypes_([&]() {
        PhenotypeList phenotypes;
        phenotypes.reserve(n_pheno);
        for (std::size_t pheno = 0; pheno < n_pheno; pheno++) {
          phenotypes.emplace_back(
              buf_,
              arch_,
              pheno_names[pheno],
              config.v_h2_gen[pheno],
              config.v_h2_env[pheno],
              config.v_h2_vert[pheno]);
        }
        return phenotypes;
      }()),
      model_(
          phenotypes_,
          config.mate_cor,
          config.n_itr,
          n_ind / 2,
          rng_,
          config.temp_init,
          config.temp_decay),
      ctx_(genome_, arch_, buf_, phenotypes_, model_),
      metrics_([&]() {
        std::vector<Metric> metrics;
        metrics.reserve(config.specs.size());
        for (const MetricSpec& spec : config.specs)
          metrics.push_back(spec.setup(ctx_));
        return metrics;
      }()) {
  if (!std::filesystem::exists(out_dir)) {
    std::filesystem::create_directory(out_dir);
  }
}

Simulation::Simulation(
    std::size_t n_gen,
    std::size_t n_ind,
    std::filesystem::path out_dir,
    std::uint64_t rng_seed,
    std::size_t n_loc,
    std::vector<double> v_maf,
    std::vector<double> v_rec,
    std::vector<double> v_mut,
    std::size_t n_pheno,
    std::vector<std::string> v_name,
    std::vector<std::size_t> v_n_loc,
    std::vector<double> v_h2_gen,
    std::vector<double> v_h2_env,
    std::vector<double> v_h2_vert,
    std::vector<double> gen_cor,
    std::vector<double> env_cor,
    std::vector<double> mate_cor,
    std::size_t n_itr,
    double tmp_init,
    double tmp_decay,
    std::vector<MetricSpec> specs,
    bool require_latent)
    : n_gen(n_gen),
      n_ind(n_ind),
      n_loc(n_loc),
      n_pheno(n_pheno),
      out_dir(out_dir),
      rng_(rng::seed_xoshiro(rng::auto_seed(rng_seed))),
      genome_(n_ind, n_loc, v_mut, v_rec, v_maf, rng_),
      arch_(
          n_pheno, n_loc, v_n_loc, v_h2_gen, v_h2_env, gen_cor, env_cor, rng_),
      buf_(n_ind, n_pheno, require_latent),
      phenotypes_([&]() {
        PhenotypeList phenotypes;
        phenotypes.reserve(n_pheno);
        for (std::size_t pheno = 0; pheno < n_pheno; pheno++) {
          phenotypes.emplace_back(
              buf_,
              arch_,
              v_name[pheno],
              v_h2_gen[pheno],
              v_h2_env[pheno],
              v_h2_vert[pheno]);
        }
        return phenotypes;
      }()),
      model_(
          phenotypes_, mate_cor, n_itr, n_ind / 2, rng_, tmp_init, tmp_decay),
      ctx_(genome_, arch_, buf_, phenotypes_, model_),
      metrics_([&]() {
        std::vector<Metric> metrics;
        metrics.reserve(specs.size());
        for (const MetricSpec& spec : specs)
          metrics.push_back(spec.setup(ctx_));
        return metrics;
      }()) {
  if (!std::filesystem::exists(out_dir))
    std::filesystem::create_directory(out_dir);
  std::cerr << "done!\n";
};

void Simulation::stream_(std::size_t gen) {
  if (gen == 0) {
    streams_.reserve(metrics_.size());

    for (std::size_t metric = 0; metric < metrics_.size(); ++metric) {
      auto stream = std::make_unique<std::ofstream>(
          out_dir / (metrics_[metric].name + ".tsv"));
      if (!stream->is_open()) {
        throw std::runtime_error(
            "could not open metric stream: " + metrics_[metric].name);
        std::cout << out_dir / (metrics_[metric].name + ".tsv") << "\n";
      }
      *stream << metrics_[metric].header() << "\n";
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
  std::cerr << "generating initial haplotypes\n";
  std::vector<std::size_t> sib_matching(ctx_.n_ind / 2);
  std::iota(sib_matching.begin(), sib_matching.end(), ctx_.n_ind / 2);

  genome_.generate_haplotypes();
  genome_.compute_mafs();
  genome_.compute_stats();

  for (std::size_t gen = 0; gen < n_gen; gen++) {
    std::cerr << "generation " << gen << "\n";
    std::cerr << "computing mafs and stats\n";
    genome_.compute_mafs();
    genome_.compute_stats();

    std::cerr << "generating env\n";
    arch_.gen_env(buf_(ComponentType::ENVIRONMENTAL), ctx_.n_ind);

    std::cerr << "scoring phenotypes\n";
    for (Phenotype& pheno : phenotypes_) {
      pheno.score(genome_);
      pheno.compute_stats();
      if (gen == 0) pheno.transmit_vert(sib_matching);
    }

    if (buf_.has_lat()) buf_.score_latent(model_.cor_U, model_.cor_VT);

    std::cerr << "mating\n";
    model_.init_state();
    model_.update(phenotypes_);
    std::vector<std::size_t> opt_matching = model_.match();

    for (Phenotype& pheno : phenotypes_) pheno.transmit_vert(opt_matching);

    std::cerr << "streaming\n";
    stream_(gen);

    std::cerr << "updating\n";
    genome_.transpose();
    genome_.update(opt_matching);
    genome_.transpose();
  }
}

void run_simulations(
    const SimulationConfig& config,
    std::size_t n_replicates,
    std::size_t n_threads) {
  // ensure the base directory exists
  if (!std::filesystem::exists(config.out_dir)) {
    std::filesystem::create_directory(config.out_dir);
  }

  // run multithreaded replicate simulations
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> pool;
  pool.reserve(n_threads);

  for (std::size_t ii = 0; ii < n_threads; ++ii) {
    pool.emplace_back([&]() {
      while (true) {
        std::size_t rep_id = next.fetch_add(1);
        if (rep_id >= n_replicates) return;

        // adjust the output directory path for the replicate
        std::filesystem::path rep_dir =
            config.out_dir / std::format("rep_{:03}", rep_id);

        // adjust replicate rng seed  
        const std::uint64_t PHI = 0x9E3779B97F4A7C15ull;
        std::uint64_t rep_seed = config.rng_seed + PHI * rep_id;

        Simulation rep(config, rep_dir, rep_seed);

        rep.run();
      }
    });
  }

  for (std::thread& thread : pool) thread.join();
}

}  // namespace amsim
