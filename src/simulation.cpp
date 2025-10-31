#include <amsim/metric.h>
#include <amsim/simulation.h>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace amsim {

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
    bool require_latent,
    std::optional<std::uint64_t> rep_id)
    : n_gen(n_gen),
      out_dir([&]() {
        if (rep_id) out_dir = out_dir / (std::format("rep_{:03d}", *rep_id));
        return out_dir;
      }()),
      status_(SimulationStatus::READY),
      rng_([&]() {
        if (rep_id) {
          const std::uint64_t gr_mult = 0x9E3779B97F4A7C15ULL;
          rng_seed = rng_seed + gr_mult * (*rep_id);
        }
        return rng::seed_xoshiro(rng::auto_seed(rng_seed));
      }()),
      genome_(n_ind, n_loc, v_mut, v_rec, v_maf, rng_),
      arch_([&]() {
        PhenoArch arch(
            n_pheno,
            n_loc,
            v_n_loc,
            v_h2_gen,
            v_h2_env,
            gen_cor,
            env_cor,
            rng_);
        arch.optim_arch(1e4);
        return arch;
      }()),
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

  ctx_.genome.generate_haplotypes();
  ctx_.genome.compute_mafs();
  ctx_.genome.compute_stats();

  for (std::size_t gen = 0; gen < n_gen; gen++) {
    std::cerr << "generation " << gen << "\n";
    std::cerr << "computing mafs and stats\n";
    ctx_.genome.compute_mafs();
    ctx_.genome.compute_stats();

    std::cerr << "generating env\n";
    ctx_.arch.gen_env(ctx_.buf(ComponentType::ENVIRONMENTAL), ctx_.n_ind);

    std::cerr << "scoring phenotypes\n";
    for (Phenotype& pheno : ctx_.phenotypes) {
      pheno.score(ctx_.genome);
      pheno.compute_stats();
      if (gen == 0) pheno.transmit_vert(sib_matching);
    }

    if (ctx_.buf.has_lat())
      ctx_.buf.score_latent(ctx_.model.cor_U, ctx_.model.cor_VT);

    std::cerr << "mating\n";
    ctx_.model.init_state();
    ctx_.model.update(ctx_.phenotypes);
    std::vector<std::size_t> opt_matching = ctx_.model.match();

    for (Phenotype& pheno : ctx_.phenotypes) pheno.transmit_vert(opt_matching);

    std::cerr << "streaming\n";
    stream_(gen);

    std::cerr << "updating\n";
    ctx_.genome.transpose();
    ctx_.genome.update(opt_matching);
    ctx_.genome.transpose();
  }
}

}  // namespace amsim
