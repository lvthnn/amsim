#ifndef AMSIMCPP_METRIC_H
#define AMSIMCPP_METRIC_H

#pragma once
#include <string>
#include <fstream>
#include <utility>

#include <amsim/genome.h>
#include <amsim/phenotype.h>

namespace amsim {
  using MetricFunc = std::function<std::vector<double>(std::vector<Phenotype>&, Genome&)>;

  class Metric {
  public:
    Metric(MetricFunc f, const std::string& file, std::string name,
           const std::size_t n_rows, const std::size_t n_cols = 1);
    void stream(std::vector<Phenotype>& phenotypes, Genome& genome);

  private:
    MetricFunc f_;
    std::size_t it_;
    std::ofstream file_;
    const std::string name_;
    const std::size_t n_rows_;
    const std::size_t n_cols_;
    std::vector<double> buf_;
  };

  namespace metrics {
    Metric make_metric(MetricFunc f, const std::string& name,
                       const std::string& file, const std::size_t n_rows,
                       const std::size_t n_cols = 1);

    namespace genome {
      inline std::vector<double> floc_maf(std::vector<Phenotype> &phenotypes, Genome &genome) { return genome.v_lmaf(); }
      inline std::vector<double> floc_mean(std::vector<Phenotype> &phenotypes, Genome &genome) { return genome.v_lmean(); }
      inline std::vector<double> floc_var(std::vector<Phenotype> &phenotype, Genome &genome) { return genome.v_lvar(); }
    }

    namespace phenome {
      MetricFunc fcomp_cor(ComponentType type);
    }

    Metric comp_cor(const std::string& file, const std::size_t n_pheno, ComponentType type);
    inline auto loc_maf(const std::string& file, const std::size_t n_loc) { return make_metric(genome::floc_maf, "loc_mafs", file, n_loc); }
    inline auto loc_mean(const std::string& file, const std::size_t n_loc) { return make_metric(genome::floc_mean, "loc_means", file, n_loc); }
    inline auto loc_var(const std::string& file, const std::size_t n_loc) { return make_metric(genome::floc_var, "loc_vars", file, n_loc); }
  }
}

#endif //AMSIMCPP_METRIC_H
