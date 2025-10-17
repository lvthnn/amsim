#include <fstream>
#include <filesystem>

#include <amsim/simulation.h>

namespace amsim {

  Simulation::Simulation(std::size_t n_gen, std::string out_dir, Genome &genome,
                         PhenoArch &arch, PhenoBuf &buf, AssortativeModel &model,
                         std::vector<Metric> &metrics)
    : n_gen_(n_gen),
      out_dir_(std::move(out_dir)),
      genome_(genome),
      arch_(arch),
      buf_(buf),
      model_(model),
      metrics_(std::move(metrics)) {
    if (!std::filesystem::exists(out_dir_))
      std::filesystem::create_directory(out_dir_);

    for (std::size_t metric = 0; metric < metrics_.size(); metric++) {
      std::ofstream file_metric; 
    }
  }

  void Simulation::stream_() {

  }

  void Simulation::run() {
    genome_.generate_haplotypes();

    for (std::size_t gen = 0; gen < n_gen_; gen++) {
      for (const auto pheno : phenotypes_)
        pheno.get().score(genome_);

      model_.update(phenotypes_);

      std::vector<std::size_t> opt_matching = model_.match();

      // for (const auto metric : metrics_) {
      //   metric.get().stream(phenotypes_, genome_);
      // }

      genome_.update(opt_matching);
    }
  }
}
