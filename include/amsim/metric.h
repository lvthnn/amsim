#ifndef AMSIMCPP_METRIC_H
#define AMSIMCPP_METRIC_H

#pragma once

#include <amsim/simulation_context.h>

#include <string>
#include <vector>

namespace amsim {

class Metric {
 public:
  Metric(
      MetricFunc f,
      const std::string& name,
      const std::size_t n_rows,
      const std::size_t n_cols = 1,
      std::vector<std::string> labels = {},
      const bool require_lat = false);

  const std::string name;
  const std::size_t n_rows;
  const std::size_t n_cols;
  const bool require_lat;

  std::string header();
  std::string stream(const SimulationContext& ctx);

 private:
  MetricFunc f_;
  std::vector<double> buf_;
  std::vector<std::string> labels_;
};

namespace metrics {

namespace genome {
inline std::vector<double> floc_mean(const SimulationContext& ctx) {
  return ctx.genome.v_lmean();
}

inline std::vector<double> floc_var(const SimulationContext& ctx) {
  return ctx.genome.v_lvar();
}

inline std::vector<double> floc_maf(const SimulationContext& ctx) {
  return ctx.genome.v_lmaf();
}
}  // namespace genome

namespace phenome {
std::vector<double> f_pheno_h2(const SimulationContext&);
std::vector<double> f_latent_h2(const SimulationContext&);
MetricFunc f_comp_cor(ComponentType type);
MetricFunc f_comp_xcor(ComponentType type);
MetricFunc f_comp_mean(ComponentType type);
MetricFunc f_comp_var(ComponentType type);
MetricFunc f_latent_comp_cor(ComponentType type);
MetricFunc f_latent_comp_xcor(ComponentType type);
MetricFunc f_latent_comp_mean(ComponentType type);
MetricFunc f_latent_comp_var(ComponentType type);
}  // namespace phenome

}  // namespace metrics

}  // namespace amsim

#endif  // AMSIMCPP_METRIC_H
