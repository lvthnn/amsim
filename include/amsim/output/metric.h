#pragma once

#include <amsim/params.h>
#include <amsim/state.h>

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace amsim {

class Metric {
 public:
  Metric(
      MetricFunc f,
      std::string name,
      std::size_t n_rows,
      std::size_t n_cols = 1,
      std::vector<std::string> labels = {});

  const std::string name;
  const std::size_t n_rows;
  const std::size_t n_cols;

  std::string header();

  std::string stream(const State& s, const Params& p);

 private:
  MetricFunc f_;                     ///< Metric computation function
  std::vector<double> data_;         ///< Data buffer
  std::vector<std::string> labels_;  ///< Column labels
};

/// Namespace for built-in metric functions
namespace metrics {

namespace genome {

Eigen::VectorXd floc_mean(const State& s, const Params& /*p*/) {
  return s.geno.v_lmean();
}

Eigen::VectorXd floc_var(const State& s, const Params& /*p*/) {
  return s.geno.v_lvar();
}

Eigen::VectorXd floc_maf(const State& s, const Params& /*p*/) {
  return s.geno.v_lmaf();
}

}  // namespace genome

namespace phenome {

std::vector<double> f_pheno_h2(const State& s, const Params& p);

MetricFunc f_comp_cor(amsim::phenome::ComponentType type);

MetricFunc f_comp_cor(
    amsim::phenome::ComponentType type_l, amsim::phenome::ComponentType type_r);

MetricFunc f_comp_xcor(amsim::phenome::ComponentType type);

MetricFunc f_comp_xcor(
    amsim::phenome::ComponentType type_l, amsim::phenome::ComponentType type_r);

MetricFunc f_comp_mean(amsim::phenome::ComponentType type);

MetricFunc f_comp_var(amsim::phenome::ComponentType type);

}  // namespace phenome

}  // namespace metrics

}  // namespace amsim
