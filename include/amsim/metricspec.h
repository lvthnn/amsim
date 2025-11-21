#ifndef AMSIMCPP_METRICSPEC_H
#define AMSIMCPP_METRICSPEC_H

#include <amsim/metric.h>
#include <amsim/simulation_context.h>

#include <functional>
#include <optional>
#include <string>

namespace amsim {

using MetricSetup = std::function<Metric(const SimulationContext& ctx)>;

class MetricSpec {
 public:
  MetricSpec(
      std::string name,
      MetricFunc func,
      MetricSetup setup,
      std::optional<bool> require_lat_ = false)
      : require_lat(require_lat_),
        name_(std::move(name)),
        func_(std::move(func)),
        setup_(std::move(setup)) {}

  Metric setup(SimulationContext& ctx) const { return setup_(ctx); };

  const bool require_lat;

 private:
  std::string name_;
  MetricFunc func_;
  MetricSetup setup_;
};

MetricSpec pheno_h2();
MetricSpec pheno_comp_cor(ComponentType type);
MetricSpec pheno_comp_xcor(ComponentType type);
MetricSpec pheno_comp_mean(ComponentType type);
MetricSpec pheno_comp_var(ComponentType type);
MetricSpec pheno_latent_h2();
MetricSpec pheno_latent_comp_cor(ComponentType type);
MetricSpec pheno_latent_comp_xcor(ComponentType type);
MetricSpec pheno_latent_comp_mean(ComponentType type);
MetricSpec pheno_latent_comp_var(ComponentType type);

}  // namespace amsim

#endif  // AMSIMCPP_METRICSPEC_H
