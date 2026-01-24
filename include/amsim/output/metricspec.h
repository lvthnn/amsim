#pragma once

#include <amsim/data/phenome.h>

#include <amsim/state.h>
#include <amsim/params.h>

#include <amsim/output/metric.h>

#include <functional>
#include <optional>
#include <string>

namespace amsim {

/// @brief Function type for setting up a metric from context
using MetricSetup = std::function<Metric(const State& s, const Params& p)>;
using MetricFunc = std::function<std::vector<double>(const State& s, const Params& p)>;

/// @brief Specification for creating a metric
///
/// MetricSpec stores the information needed to create a Metric instance,
/// including the computation function, setup function, and metadata.
class MetricSpec {
 public:
  /// @brief Construct a MetricSpec
  ///
  /// @param name Metric name
  /// @param func Function computing metric values
  /// @param setup Function to create Metric from context
  /// @param require_lat_ Whether metric requires latent phenotypes
  MetricSpec(
      std::string name,
      MetricFunc func,
      MetricSetup setup,
      std::optional<bool> require_lat_ = false)
      : require_lat(require_lat_),
        name_(std::move(name)),
        func_(std::move(func)),
        setup_(std::move(setup)) {}

  /// @brief Create a Metric instance from context
  /// @param ctx Simulation context
  /// @return Configured Metric instance
  Metric setup(const State& s, const Params& p) const { return setup_(s, p); };

  const bool require_lat;  ///< Whether latent phenotypes are required

 private:
  std::string name_;   ///< Metric name
  MetricFunc func_;    ///< Metric computation function
  MetricSetup setup_;  ///< Metric setup function
};

// @brief Create spec for locus mean metric
MetricSpec loc_mean();

// @brief Create spec for locus variance metric
MetricSpec loc_var();

// @brief Create spec for locus MAF metric
MetricSpec loc_maf();

/// @brief Create spec for phenotype heritability metric
MetricSpec pheno_h2();

/// @brief Create spec for phenotype component correlation metric
/// @param type Component type
MetricSpec pheno_comp_cor(ComponentType type);

/// @brief Create spec for phenotype cross-component correlation metric
/// @param type Left component type
/// @param type Right component type
MetricSpec pheno_comp_cor(ComponentType type_l, ComponentType type_r);

/// @brief Create spec for phenotype component cross-correlation metric
/// @param type Component type
MetricSpec pheno_comp_xcor(ComponentType type);

/// @brief Create spec for phenotype cross-component cross-correlation metric
/// @param type Left component type
/// @param type Right component type
MetricSpec pheno_comp_xcor(ComponentType type_l, ComponentType type_r);

/// @brief Create spec for phenotype component mean metric
/// @param type Component type
MetricSpec pheno_comp_mean(ComponentType type);

/// @brief Create spec for phenotype component variance metric
/// @param type Component type
MetricSpec pheno_comp_var(ComponentType type);

/// @brief Create spec for latent phenotype heritability metric
MetricSpec pheno_latent_h2();

/// @brief Create spec for latent phenotype component correlation metric
/// @param type Component type
MetricSpec pheno_latent_comp_cor(ComponentType type);

/// @brief Create spec for latent phenotype component cross-correlation metric
/// @param type Component type
MetricSpec pheno_latent_comp_xcor(ComponentType type);

/// @brief Create spec for latent phenotype component mean metric
/// @param type Component type
MetricSpec pheno_latent_comp_mean(ComponentType type);

/// @brief Create spec for latent phenotype component variance metric
/// @param type Component type
MetricSpec pheno_latent_comp_var(ComponentType type);

}  // namespace amsim
