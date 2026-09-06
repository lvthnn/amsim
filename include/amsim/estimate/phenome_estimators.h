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

#include <amsim/core/params.h>
#include <amsim/core/state.h>
#include <amsim/core/utils.h>
#include <amsim/estimate/estimator.h>
#include <amsim/estimate/strategy.h>

#include <Eigen/Dense>
#include <optional>

namespace amsim {

namespace details {

class EstimatorHeritability : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorHeritability(const Params& params)
      : PopulationEstimatorStrategy(
            "pheno_h2", {}, params.pheno.names, params.pheno.n_pheno) {
    n_pheno_ = params.pheno.n_pheno;
  }

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno) = state.pheno().comp_var(pheno, Component::Genetic) /
                     state.pheno().comp_var(pheno, Component::Total);
  }

 private:
  std::size_t n_pheno_;
};

class EstimatorComponentMean : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentMean(const Params& params, Component type)
      : PopulationEstimatorStrategy(
            "pheno_" + to_string(type) + "_mean",
            {},
            params.pheno.names,
            params.pheno.n_pheno),
        type_(type),
        n_pheno_(params.pheno.n_pheno) {}

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno) = state.pheno().comp_mean(pheno, type_);
  }

 private:
  Component type_;
  std::size_t n_pheno_;
};

class EstimatorComponentVar : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentVar(const Params& params, Component type)
      : PopulationEstimatorStrategy(
            "pheno_" + to_string(type) + "_var",
            {},
            params.pheno.names,
            params.pheno.n_pheno),
        type_(type),
        n_pheno_(params.pheno.n_pheno) {}

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno) = state.pheno().comp_var(pheno, type_);
  }

 private:
  Component type_;
  std::size_t n_pheno_;
};

class EstimatorComponentCor : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentCor(
      const Params& params, Component type_l, std::optional<Component> type_r)
      : PopulationEstimatorStrategy(
            "pheno_" + to_string(type_l) + "_" +
                to_string(type_r.value_or(type_l)) + "_cor",
            utils::vector_suffix(params.pheno.names, "_" + to_string(type_l)),
            utils::vector_suffix(
                params.pheno.names, "_" + to_string(type_r.value_or(type_l))),
            params.pheno.n_pheno,
            params.pheno.n_pheno),
        n_ind_(params.global.n_ind),
        n_pheno_(params.pheno.n_pheno),
        type_l_(type_l),
        type_r_(type_r.value_or(type_l)),
        std_l_(n_ind_, n_pheno_),
        std_r_(n_ind_, n_pheno_) {}

  void compute(const State& state) override {
    std_l_ = utils::standardise(state.pheno()(type_l_));
    std_r_ = utils::standardise(state.pheno()(type_r_));
    data_ = (std_l_.transpose() * std_r_) / static_cast<double>(n_ind_);
  }

 private:
  std::size_t n_ind_;
  std::size_t n_pheno_;
  Component type_l_;
  Component type_r_;
  Eigen::MatrixXd std_l_;
  Eigen::MatrixXd std_r_;
};

class EstimatorComponentCov : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentCov(
      const Params& params, Component type_l, std::optional<Component> type_r)
      : PopulationEstimatorStrategy(
            "pheno_" + to_string(type_l) + "_" +
                to_string(type_r.value_or(type_l)) + "_cov",
            utils::vector_suffix(params.pheno.names, "_" + to_string(type_l)),
            utils::vector_suffix(
                params.pheno.names, "_" + to_string(type_r.value_or(type_l))),
            params.pheno.n_pheno,
            params.pheno.n_pheno),
        n_ind_(params.global.n_ind),
        n_pheno_(params.pheno.n_pheno),
        type_l_(type_l),
        type_r_(type_r.value_or(type_l)),
        centred_l_(n_ind_, n_pheno_),
        centred_r_(n_ind_, n_pheno_) {}

  void compute(const State& state) override {
    centred_l_ = utils::standardise(state.pheno()(type_l_), true, true, false);
    centred_r_ = utils::standardise(state.pheno()(type_r_), true, true, false);
    data_ =
        (centred_l_.transpose() * centred_r_) / static_cast<double>(n_ind_ - 1);
  }

 private:
  std::size_t n_ind_;
  std::size_t n_pheno_;
  Component type_l_;
  Component type_r_;
  Eigen::MatrixXd centred_l_;
  Eigen::MatrixXd centred_r_;
};

}  // namespace details

inline PopulationEstimator PopulationHeritability() {
  return PopulationEstimator{
      .name = "heritability", .fn = [](const Params& params) {
        return std::make_unique<details::EstimatorHeritability>(params);
      }};
}

inline PopulationEstimator PopulationComponentMean(
    Component type = Component::Total) {
  return PopulationEstimator{
      .name = "pheno-mean-" + to_string(type),
      .fn = [type](const Params& params) {
        return std::make_unique<details::EstimatorComponentMean>(params, type);
      }};
}

inline PopulationEstimator PopulationComponentVar(
    Component type = Component::Total) {
  return PopulationEstimator{
      .name = "pheno-var-" + to_string(type),
      .fn = [type](const Params& params) {
        return std::make_unique<details::EstimatorComponentVar>(params, type);
      }};
}

inline PopulationEstimator PopulationComponentCor(
    Component type_l = Component::Total,
    std::optional<Component> type_r = std::nullopt) {
  return PopulationEstimator{
      .name = "pheno-cor-" + to_string(type_l) + "-" +
              to_string(type_r.value_or(type_l)),
      .fn = [type_l, type_r](const Params& params) {
        return std::make_unique<details::EstimatorComponentCor>(
            params, type_l, type_r);
      }};
}

inline PopulationEstimator PopulationComponentCov(
    Component type_l = Component::Total,
    std::optional<Component> type_r = std::nullopt) {
  return PopulationEstimator{
      .name = "pheno-cov-" + to_string(type_l) + "-" +
              to_string(type_r.value_or(type_l)),
      .fn = [type_l, type_r](const Params& params) {
        return std::make_unique<details::EstimatorComponentCov>(
            params, type_l, type_r);
      }};
}

}  // namespace amsim
