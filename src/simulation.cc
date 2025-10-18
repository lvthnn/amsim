#include <filesystem>
#include <stdexcept>
#include <iostream>

#include <amsim/simulation.h>

namespace amsim {
  Simulation::Simulation(// SIMULATION PARAMETERS
                          std::size_t n_gen,
                          std::size_t n_ind,
                          std::filesystem::path out_dir,
                          std::uint64_t rng_seed,
                          // GENOME PARAMETERS
                          std::size_t n_loc,
                          std::vector<double> v_maf,
                          std::vector<double> v_rec,
                          std::vector<double> v_mut,
                          // PHENOME PARAMETERS
                          std::size_t n_pheno,
                          std::vector<std::string> v_name,
                          std::vector<std::size_t> v_n_loc,
                          std::vector<double> v_h2_gen,
                          std::vector<double> v_h2_env,
                          std::vector<double> v_h2_vert,
                          std::vector<double> gen_cor,
                          std::vector<double> env_cor,
                          // MATING MODEL PARAMETERS
                          std::vector<double> mate_cor,
                          std::size_t n_itr,
                          double tmp_init,
                          double tmp_decay,
                          std::vector<Metric> metrics)
		: n_gen(n_gen),
			out_dir(out_dir),
			status_(SimulationStatus::READY),
      rng_(rng::seed_xoshiro(rng::auto_seed(rng_seed))),
			genome_(n_ind, n_loc, v_mut, v_rec, v_maf, rng_),
      arch_([&](){
        PhenoArch arch(n_pheno, n_loc, v_n_loc, v_h2_gen, v_h2_env, 
                      gen_cor, env_cor, rng_);
        arch.optim_arch(1e4);
        return arch;
      }()),
      buf_(n_ind, n_pheno),
      phenotypes_([&](){
        PhenotypeList phenotypes; 
        phenotypes.reserve(n_pheno);
        for (std::size_t pheno = 0; pheno < n_pheno; pheno++) {
          phenotypes.emplace_back(buf_, arch_, v_name[pheno], v_h2_gen[pheno],
                                  v_h2_env[pheno], v_h2_vert[pheno]);
        }
        return phenotypes;
      }()),
      model_(phenotypes_, mate_cor, n_itr, n_ind / 2, rng_, tmp_init, tmp_decay),
		  metrics_(std::move(metrics)),
      ctx(genome_, arch_, buf_, phenotypes_, model_) {
    streams_.resize(metrics_.size());
    for (std::size_t metric = 0; metric < metrics_.size(); metric++) {
      auto stream = std::make_unique<std::ofstream>(out_dir / (metrics_[metric].name + ".tsv"));
      if (!stream->is_open())
        throw std::runtime_error("could not open metric stream");
      streams_.push_back(std::move(stream));
    }
  };

  void Simulation::stream_(std::size_t gen) {
    for (std::size_t metric = 0; metric < metrics_.size(); metric++) {
      if (!streams_[metric]) streams_.resize(metrics_.size());
      if (!streams_[metric] || !streams_[metric]->is_open()) {
        streams_[metric] = std::make_unique<std::ofstream>(
          out_dir / (metrics_[metric].name + ".tsv")
        );
        if (!streams_[metric]->is_open()) {
          throw std::runtime_error("Could not open metric file: " + metrics_[metric].name);
        }
        *streams_[metric] << metrics_[metric].header() << "\n";
      }
      *streams_[metric] << std::to_string(gen + 1) << "\t";
      *streams_[metric] << metrics_[metric].stream(ctx) << "\n";
    }
  }


  void Simulation::run() {
    ctx.genome.generate_haplotypes();
    ctx.genome.compute_mafs();
    ctx.genome.compute_stats();

    for (std::size_t gen = 0; gen < n_gen; gen++) {
      ctx.genome.compute_mafs();
      ctx.genome.compute_stats();

      ctx.arch.gen_env(ctx.buf(0, ComponentType::ENVIRONMENTAL), ctx.buf.n_ind());

      for (Phenotype& pheno : ctx.phenotypes) {
        pheno.score(ctx.genome);
        pheno.compute_stats();
      } 

      ctx.model.update(ctx.phenotypes);
      std::vector<std::size_t> opt_matching = ctx.model.match();
       ctx.model.display_cor();

      stream_(gen);

      ctx.genome.transpose();

      uint64_t before_H0 = ctx.genome.H0()(0, 0);
      uint64_t before_H1 = ctx.genome.H1()(0, 0);

      ctx.genome.update(opt_matching);

      uint64_t after_H0 = ctx.genome.H0()(0, 0);
      uint64_t after_H1 = ctx.genome.H1()(0, 0);

      if (before_H0 == after_H0 && before_H1 == after_H1) {
        std::cerr << "WARNING: Haplotypes unchanged after update!\n";
      } else {
        std::cerr << "Haplotypes changed: H0 " << before_H0 << " -> " << after_H0 << "\n";
      }

      ctx.genome.transpose();
    }
    std::cerr << "All generations complete\n";

  }
}
