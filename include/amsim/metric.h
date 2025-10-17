#ifndef AMSIMCPP_METRIC_H
#define AMSIMCPP_METRIC_H

#pragma once
#include <string>
#include <optional>

#include <amsim/genome.h>
#include <amsim/phenotype.h>

namespace amsim {
  using MetricFunc = std::function<std::vector<double>(PhenotypeList&, Genome&)>;

  class Metric {
  public:
    Metric(MetricFunc f, const std::string &name,
           const std::size_t n_rows, const std::size_t n_cols = 1,
           std::optional<std::vector<std::string>> el_names = std::nullopt);

    inline const std::string name() const noexcept { return name_; }
    inline const std::size_t n_rows() const noexcept { return n_rows_; }
    inline const std::size_t n_cols() const noexcept { return n_cols_; }
    inline const bool named() const noexcept { return named_; }
    inline const std::vector<std::string> names() const noexcept { return el_names_; }
    void stream(PhenotypeList &phenotypes, Genome &genome);

  private:
    MetricFunc f_;
    const std::string name_;
    const std::size_t n_rows_;
    const std::size_t n_cols_;
    std::vector<double> buf_;

    bool named_ = false;
    std::vector<std::string> el_names_;

    std::size_t it_ = 0;
  };

  namespace metrics {
    Metric make_metric(MetricFunc f, const std::string& name,
                       const std::size_t n_rows,
                       const std::size_t n_cols = 1,
                       std::optional<std::vector<std::string>> el_names = std::nullopt);

    namespace genome {
      inline std::vector<double> floc_mean(PhenotypeList &phenotypes, Genome &genome) { return genome.v_lmean(); }
      inline std::vector<double> floc_var(PhenotypeList &phenotypes, Genome &genome) { return genome.v_lvar(); }
      inline std::vector<double> floc_maf(PhenotypeList &phenotypes, Genome &genome) { return genome.v_lmaf(); }
    }

    namespace phenome {
      MetricFunc fcomp_cor(ComponentType type);
    }

    Metric pheno_h2(const std::size_t n_pheno);
    Metric comp_cor(const std::size_t n_pheno, ComponentType type);
    inline Metric loc_maf(const std::size_t n_loc) { return make_metric(genome::floc_maf, "loc_mafs", n_loc); }
    inline Metric loc_mean(const std::size_t n_loc) { return make_metric(genome::floc_mean, "loc_means", n_loc); }
    inline Metric loc_var(const std::size_t n_loc) { return make_metric(genome::floc_var, "loc_vars", n_loc); }
  }

  using MetricList = std::vector<std::reference_wrapper<Metric>>;
}

#endif //AMSIMCPP_METRIC_H
