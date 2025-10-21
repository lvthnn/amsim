#ifndef AMSIMCPP_METRIC_H
#define AMSIMCPP_METRIC_H

#pragma once
#include <string>
#include <optional>
#include <vector>

#include <amsim/simulation_context.h>

namespace amsim {
  class Metric {
  public:
    Metric(MetricFunc f, const std::string &name,
           const std::size_t n_rows, const std::size_t n_cols = 1,
           const bool require_lat = false,
           std::optional<std::vector<std::string>> el_names = std::nullopt);

    const std::string name;
    const std::size_t n_rows;
    const std::size_t n_cols;
    const bool named;
    const std::vector<std::string> names;
    const bool require_lat;

    std::string header();
    std::string stream(const SimulationContext &ctx);

  private:
    MetricFunc f_;
    std::vector<double> buf_;
  };

  namespace metrics {
    Metric make_metric(MetricFunc f, const std::string& name,
                       const std::size_t n_rows,
                       const std::size_t n_cols = 1,
                       const bool require_lat = false,
                       std::optional<std::vector<std::string>> el_names = std::nullopt);

    namespace genome {
      inline std::vector<double> floc_mean(const SimulationContext &ctx) {
        return ctx.genome.v_lmean();
      }

      inline std::vector<double> floc_var(const SimulationContext &ctx) {
        return ctx.genome.v_lvar();
      }

      inline std::vector<double> floc_maf(const SimulationContext &ctx) {
        return ctx.genome.v_lmaf();
      }
    }

    namespace phenome {
      MetricFunc f_comp_cor(ComponentType type);
    }

    Metric pheno_h2(const std::size_t n_pheno);
    Metric comp_mean(const std::size_t n_pheno, ComponentType type);
    Metric comp_var(const std::size_t n_pheno, ComponentType type);
    Metric comp_cor(const std::size_t n_pheno, ComponentType type);
    Metric comp_xcor(const std::size_t n_pheno, ComponentType type);
    Metric latent_h2(const std::size_t n_pheno);
    Metric latent_comp_cor(const std::size_t n_pheno, ComponentType type);
    Metric latent_comp_xcor(const std::size_t n_pheno, ComponentType type);
    Metric latent_comp_var(const std::size_t n_pheno, ComponentType type);

    inline Metric loc_maf(const std::size_t n_loc) {
      return make_metric(genome::floc_maf, "loc_mafs", n_loc);
    }

    inline Metric loc_mean(const std::size_t n_loc) {
      return make_metric(genome::floc_mean, "loc_means", n_loc);
    }

    inline Metric loc_var(const std::size_t n_loc) {
      return make_metric(genome::floc_var, "loc_vars", n_loc);
    }
  }
}

#endif //AMSIMCPP_METRIC_H
