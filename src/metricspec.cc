#include <amsim/component_type.h>
#include <amsim/metric.h>
#include <amsim/metricspec.h>
#include <amsim/simulation_context.h>

#include <string>
#include <vector>

namespace amsim {

namespace labels {

std::vector<std::string> label_matrix(
    const std::vector<std::string>& names,
    const std::optional<std::vector<std::string>>& suffix = std::nullopt) {
  std::size_t n = names.size();
  std::vector<std::string> labels(n * n);
  for (std::size_t i = 0; i < names.size(); ++i)
    for (std::size_t j = 0; j < names.size(); ++j)
      labels[(i * n) + j] = (suffix) ? names[i] + "_" + (*suffix)[0] +
                                           "::" + names[j] + "_" + (*suffix)[1]
                                     : labels[(i * n) + j] =
                                           names[i] + "::" + names[j];
  return labels;
}
}  // namespace labels

MetricSpec loc_mean() {
  std::string name = "loc_mean";
  MetricFunc func = metrics::genome::floc_mean;
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_loc;
    const std::size_t n_cols = 1;

    std::vector<std::string> labels(ctx.n_loc);
    for (std::size_t el = 0; el < ctx.n_loc; ++el)
      labels[el] = "loc" + std::to_string(el + 1);

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup);
}

MetricSpec loc_var() {
  std::string name = "loc_var";
  MetricFunc func = metrics::genome::floc_var;
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_loc;
    const std::size_t n_cols = 1;

    std::vector<std::string> labels(ctx.n_loc);
    for (std::size_t el = 0; el < ctx.n_loc; ++el)
      labels[el] = "loc" + std::to_string(el + 1);

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup);
}

MetricSpec loc_maf() {
  std::string name = "loc_maf";
  MetricFunc func = metrics::genome::floc_maf;
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_loc;
    const std::size_t n_cols = 1;

    std::vector<std::string> labels(ctx.n_loc);
    for (std::size_t el = 0; el < ctx.n_loc; ++el)
      labels[el] = "loc" + std::to_string(el + 1);

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup);
}

MetricSpec pheno_h2() {
  std::string name = "pheno_h2";
  MetricFunc func = metrics::phenome::f_pheno_h2;
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = 1;
    const std::vector<std::string> labels = ctx.pheno_names;

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup);
}

MetricSpec pheno_comp_cor(ComponentType type) {
  std::string name = "pheno_" + to_string(type) + "_cor";
  MetricFunc func = metrics::phenome::f_comp_cor(type);
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = ctx.n_pheno;
    const std::vector<std::string> labels =
        labels::label_matrix(ctx.pheno_names);

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup);
}

MetricSpec pheno_comp_cor(ComponentType type_l, ComponentType type_r) {
  std::string name =
      "pheno_" + to_string(type_l) + "_" + to_string(type_r) + "_cor";
  std::vector<std::string> suffixes({to_string(type_l), to_string(type_r)});
  MetricFunc func = metrics::phenome::f_comp_cor(type_l, type_r);
  MetricSetup setup = [func, name, suffixes](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = ctx.n_pheno;
    const std::vector<std::string> labels =
        labels::label_matrix(ctx.pheno_names, suffixes);

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup);
}

MetricSpec pheno_comp_xcor(ComponentType type) {
  std::string name = "pheno_" + to_string(type) + "_xcor";
  MetricFunc func = metrics::phenome::f_comp_xcor(type);
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = ctx.n_pheno;
    const std::vector<std::string> labels = labels::label_matrix(
        ctx.pheno_names, std::vector<std::string>({"male", "female"}));

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup);
}

MetricSpec pheno_comp_mean(ComponentType type) {
  std::string name = "pheno_" + to_string(type) + "_mean";
  MetricFunc func = metrics::phenome::f_comp_mean(type);
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = 1;
    const std::vector<std::string> labels = ctx.pheno_names;

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup);
}

MetricSpec pheno_comp_var(ComponentType type) {
  std::string name = "pheno_" + to_string(type) + "_var";
  MetricFunc func = metrics::phenome::f_comp_var(type);
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = 1;
    const std::vector<std::string> labels = ctx.pheno_names;

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup);
}

MetricSpec pheno_latent_h2() {
  std::string name = "pheno_h2";
  MetricFunc func = metrics::phenome::f_pheno_h2;
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = 1;
    const std::vector<std::string> labels = ctx.pheno_names;

    return Metric(func, name, n_rows, n_cols, labels);
  };
  return MetricSpec(name, func, setup, true);
}

MetricSpec pheno_latent_comp_cor(ComponentType type) {
  std::string name = "latent_pheno_" + to_string(type) + "_cor";
  MetricFunc func = metrics::phenome::f_latent_comp_cor(type);
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = ctx.n_pheno;

    std::vector<std::string> latent_names(ctx.n_pheno);
    for (std::size_t el = 0; el < ctx.n_pheno; ++el)
      latent_names[el] = "L" + std::to_string(el + 1);

    const std::vector<std::string> labels = labels::label_matrix(latent_names);

    return Metric(func, name, n_rows, n_cols, labels, true);
  };
  return MetricSpec(name, func, setup, true);
}

MetricSpec pheno_latent_comp_xcor(ComponentType type) {
  std::string name = "latent_pheno_" + to_string(type) + "_xcor";
  MetricFunc func = metrics::phenome::f_latent_comp_xcor(type);
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = ctx.n_pheno;

    std::vector<std::string> latent_names(ctx.n_pheno);
    for (std::size_t el = 0; el < ctx.n_pheno; ++el)
      latent_names[el] = "L" + std::to_string(el + 1);

    const std::vector<std::string> labels = labels::label_matrix(
        latent_names, std::vector<std::string>({"male", "female"}));

    return Metric(func, name, n_rows, n_cols, labels, true);
  };
  return MetricSpec(name, func, setup, true);
}

MetricSpec pheno_latent_comp_mean(ComponentType type) {
  std::string name = "latent_pheno_" + to_string(type) + "_mean";
  MetricFunc func = metrics::phenome::f_latent_comp_mean(type);
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = 1;

    std::vector<std::string> labels(ctx.n_pheno);
    for (std::size_t el = 0; el < ctx.n_pheno; ++el)
      labels[el] = "L" + std::to_string(el + 1);

    return Metric(func, name, n_rows, n_cols, labels, true);
  };
  return MetricSpec(name, func, setup, true);
}

MetricSpec pheno_latent_comp_var(ComponentType type) {
  std::string name = "latent_pheno_" + to_string(type) + "_var";
  MetricFunc func = metrics::phenome::f_latent_comp_var(type);
  MetricSetup setup = [func, name](const SimulationContext& ctx) {
    const std::size_t n_rows = ctx.n_pheno;
    const std::size_t n_cols = 1;

    std::vector<std::string> labels(ctx.n_pheno);
    for (std::size_t el = 0; el < ctx.n_pheno; ++el)
      labels[el] = "L" + std::to_string(el + 1);

    return Metric(func, name, n_rows, n_cols, labels, true);
  };
  return MetricSpec(name, func, setup, true);
}
}  // namespace amsim
