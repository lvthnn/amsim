// This file is part of amsim, copyright (C) 2025-2026 Kári Hlynsson.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <amsim/data/component.h>
#include <amsim/estimate/sample_estimator.h>
#include <amsim/io/parse.h>
#include <amsim/sample/proband.h>

#include <algorithm>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace amsim {

// raw, pre-resolution sample declaration
struct SampleSpec {
  std::string name;
  std::string proband_type;
  std::size_t n_probands;
  std::optional<std::vector<std::string>> on;
  std::optional<std::vector<std::string>> of;
  std::optional<std::string> agg;
  std::optional<std::string> weight_function;
  std::vector<std::string> estimators;
};

// user-facing specification struct
template <Proband P>
struct Sample {
  using ProbandEnum = ProbandData<P>::ProbandEnum;

  // selection parameters — which data, from whom, and what to do
  std::string name;
  std::size_t n_probands;
  std::optional<std::vector<std::string>> on;
  std::optional<std::vector<Component>> on_components;
  ProbandEnum of = ProbandData<P>::ProbandDefault;
  Aggregator agg = ProbandData<P>::AggDefault;
  bool decompress_genotypes = false;

  // weighting method — probability of selecting based on proband aggregate
  WeightFunction weighting = Uniform();

  // estimators are declared here
  SampleEstimators<P> estimators;

  // defined out-of-line in sample/sampler.h, where build_sample_estimator's
  // full implementation (and everything it depends on) is available
  void attach_estimator(const SampleEstimatorSpec& spec);
};

using SampleVariant = std::variant<
    Sample<Proband::Individual>,
    Sample<Proband::Mate>,
    Sample<Proband::Family>>;

template <Proband P>
inline Sample<P> build_sample(const SampleSpec& spec) {
  Sample<P> sample;
  sample.name = spec.name;
  sample.n_probands = spec.n_probands;
  if (spec.on.has_value()) sample.on = spec.on.value();
  if (spec.of.has_value()) sample.of = parse_proband<P>(spec.of.value());
  if (spec.agg.has_value())
    sample.agg = Aggregator_from_string(spec.agg.value());
  if (spec.weight_function.has_value())
    sample.weighting = parse<WeightFunction>(spec.weight_function.value());
  return sample;
}

inline SampleVariant build_sample(const SampleSpec& spec) {
  if (spec.proband_type == "individual")
    return build_sample<Proband::Individual>(spec);
  if (spec.proband_type == "family") return build_sample<Proband::Family>(spec);
  if (spec.proband_type == "mate") return build_sample<Proband::Mate>(spec);
  throw std::runtime_error("Invalid Proband type " + spec.proband_type);
}

inline std::vector<SampleVariant> build_samples(
    const std::vector<SampleSpec>& spec,
    const std::vector<SampleEstimatorSpec>& estimator_spec) {
  for (const auto& se : estimator_spec) {
    if (se.type == "external") {
      if (!se.exec.has_value())
        throw std::runtime_error(
            "sample estimator '" + se.name +
            "': type 'external' requires --exec");
      if (!se.n_rows.has_value() || !se.n_cols.has_value())
        throw std::runtime_error(
            "sample estimator '" + se.name +
            "': type 'external' requires --n-rows and --n-cols");
    } else if (se.exec.has_value()) {
      throw std::runtime_error(
          "sample estimator '" + se.name +
          "': --exec is only valid for type 'external'");
    }
  }

  std::vector<SampleVariant> samples;
  for (const auto& sample_spec : spec) {
    if (sample_spec.estimators.empty())
      throw std::runtime_error(
          "sample '" + sample_spec.name + "' has no estimators declared");

    SampleVariant sample = build_sample(sample_spec);
    for (const auto& estimator : sample_spec.estimators) {
      auto it = std::ranges::find_if(
          estimator_spec,
          [&](const SampleEstimatorSpec& sample_estimator_spec) {
            return sample_estimator_spec.name == estimator;
          });

      if (it == estimator_spec.end())
        throw std::runtime_error(
            "Sample estimator " + estimator + " attached to sample " +
            sample_spec.name + " not found");

      std::visit([&](auto& sample) { sample.attach_estimator(*it); }, sample);
    }
    samples.push_back(std::move(sample));
  }
  return samples;
}

}  // namespace amsim
