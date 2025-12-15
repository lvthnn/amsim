#ifndef AMSIMCPP_METRIC_H
#define AMSIMCPP_METRIC_H

#include <amsim/simulation_context.h>

#include <string>
#include <vector>

namespace amsim {

/// @brief Metric computed during simulation
///
/// Metric wraps a function that computes values from the simulation context,
/// along with metadata for output formatting and labeling.
class Metric {
 public:
  /// @brief Construct a Metric
  ///
  /// @param f Function computing metric values
  /// @param name Metric name
  /// @param n_rows Number of output rows
  /// @param n_cols Number of output columns
  /// @param labels Column labels
  /// @param require_lat Whether metric requires latent phenotypes
  Metric(
      MetricFunc f,
      std::string name,
      std::size_t n_rows,
      std::size_t n_cols = 1,
      std::vector<std::string> labels = {},
      bool require_lat = false);

  const std::string name;    ///< Metric name
  const std::size_t n_rows;  ///< Number of output rows
  const std::size_t n_cols;  ///< Number of output columns
  const bool require_lat;    ///< Whether latent phenotypes are required

  /// @brief Generate CSV header string
  /// @return Header string with column labels
  std::string header();

  /// @brief Compute metric and format as CSV string
  /// @param ctx Simulation context
  /// @return CSV-formatted metric values
  std::string stream(const SimulationContext& ctx);

 private:
  MetricFunc f_;                     ///< Metric computation function
  std::vector<double> buf_;          ///< Value buffer
  std::vector<std::string> labels_;  ///< Column labels
};

/// Namespace for built-in metric functions
namespace metrics {

/// Genome-level metrics
namespace genome {
/// @brief Compute locus means
/// @param ctx Simulation context
/// @return Vector of locus means
inline std::vector<double> floc_mean(const SimulationContext& ctx) {
  return ctx.genome.v_lmean();
}

/// @brief Compute locus variances
/// @param ctx Simulation context
/// @return Vector of locus variances
inline std::vector<double> floc_var(const SimulationContext& ctx) {
  return ctx.genome.v_lvar();
}

/// @brief Compute locus minor allele frequencies
/// @param ctx Simulation context
/// @return Vector of locus MAFs
inline std::vector<double> floc_maf(const SimulationContext& ctx) {
  return ctx.genome.v_lmaf();
}
}  // namespace genome

/// Phenotype-level metrics
namespace phenome {
/// @brief Compute phenotype heritabilities
std::vector<double> f_pheno_h2(const SimulationContext&);

/// @brief Compute latent phenotype heritabilities
std::vector<double> f_latent_h2(const SimulationContext&);

/// @brief Create metric function for component correlations
MetricFunc f_comp_cor(ComponentType type);

/// @brief Create metric function for cross-component correlations
MetricFunc f_comp_cor(ComponentType type_l, ComponentType type_r);

/// @brief Create metric function for component cross-correlations
MetricFunc f_comp_xcor(ComponentType type);

/// @brief Create metric function for cross-component cross-correlations
MetricFunc f_comp_xcor(ComponentType type_l, ComponentType type_r);

/// @brief Create metric function for component means
MetricFunc f_comp_mean(ComponentType type);

/// @brief Create metric function for component variances
MetricFunc f_comp_var(ComponentType type);

/// @brief Create metric function for latent component correlations
MetricFunc f_latent_comp_cor(ComponentType type);

/// @brief Create metric function for latent component cross-correlations
MetricFunc f_latent_comp_xcor(ComponentType type);

/// @brief Create metric function for latent component means
MetricFunc f_latent_comp_mean(ComponentType type);

/// @brief Create metric function for latent component variances
MetricFunc f_latent_comp_var(ComponentType type);
}  // namespace phenome

}  // namespace metrics

}  // namespace amsim

#endif  // AMSIMCPP_METRIC_H
