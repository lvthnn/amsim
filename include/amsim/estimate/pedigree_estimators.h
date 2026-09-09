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

#include <amsim/core/log.h>
#include <amsim/core/params.h>
#include <amsim/core/state.h>
#include <amsim/core/utils.h>
#include <amsim/data/pedigree.h>
#include <amsim/estimate/estimator.h>
#include <amsim/estimate/strategy.h>

#include <Eigen/Dense>
#include <deque>
#include <limits>

namespace amsim {

namespace details {

class EstimatorCousinCov : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorCousinCov(
      const Params& params, std::size_t degree, Component type)
      : type_(type),
        degree_(degree),
        n_ind_(params.global.n_ind),
        n_pheno_(params.pheno.n_pheno),
        self_(params.global.n_ind * (1ULL << (2 * degree)), params.pheno.n_pheno),
        cousin_(
            params.global.n_ind * (1ULL << (2 * degree)), params.pheno.n_pheno),
        PopulationEstimatorStrategy(
            "cousin_" + std::to_string(degree) + "_" + to_string(type) + "_cov",
            utils::vector_prefix(params.pheno.names, "self_"),
            utils::vector_prefix(params.pheno.names, "cousin_"),
            params.pheno.n_pheno,
            params.pheno.n_pheno) {
    if (params.global.pedigree_max_depth < degree_ + 2)
      Log::warning(
          "Pedigree depth is insufficient to compute cousin covariance of "
          "degree " +
          std::to_string(degree_));
  }

  void compute(const State& state) override {
    if (state.pedigree.depth() < degree_ + 2) {
      data_.setConstant(std::numeric_limits<double>::quiet_NaN());
      return;
    }

    std::vector<std::vector<PedigreeNode>> paths =
        state.pedigree.find_cousins(degree_);

    auto buf = state.pheno()(type_);
    for (std::size_t ind = 0; ind < paths.size(); ++ind) {
      auto self = paths[ind].front().index;
      auto cousin = paths[ind].back().index;
      self_.row(ind) = buf.row(self);
      cousin_.row(ind) = buf.row(cousin);
    }

    self_ = utils::standardise(self_, true, true, false);
    cousin_ = utils::standardise(cousin_, true, true, false);

    data_ =
        (self_.transpose() * cousin_) / static_cast<double>(paths.size() - 1);
  }

 private:
  Component type_;
  std::size_t degree_;
  std::size_t n_ind_;
  std::size_t n_pheno_;
  Eigen::MatrixXd self_;
  Eigen::MatrixXd cousin_;
};

class EstimatorAncestorCov : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorAncestorCov(
      const Params& params,
      std::size_t degree = 1,
      Component type = Component::Genetic)
      : type_(type),
        degree_(degree),
        n_ind_(params.global.n_ind),
        n_pheno_(params.pheno.n_pheno),
        self_(params.global.n_ind * (1ULL << degree), params.pheno.n_pheno),
        ancestor_(params.global.n_ind * (1ULL << degree), params.pheno.n_pheno),
        PopulationEstimatorStrategy(
            "ancestor_" + std::to_string(degree) + "_" + to_string(type) +
                "_cov",
            utils::vector_prefix(params.pheno.names, "self_"),
            utils::vector_prefix(params.pheno.names, "ancestor_"),
            params.pheno.n_pheno,
            params.pheno.n_pheno) {}

  void compute(const State& state) override {
    syncPhenotypes(state);

    if (state.pedigree.depth() < degree_) {
      data_.setConstant(std::numeric_limits<double>::quiet_NaN());
      return;
    }

    std::vector<std::vector<PedigreeNode>> paths =
        state.pedigree.find_ancestors(degree_);

    auto self_buf = state.pheno()(type_);
    auto ancestor_buf = history_[degree_](type_);

    for (std::size_t ind = 0; ind < paths.size(); ++ind) {
      auto self = paths[ind].front().index;
      auto ancestor = paths[ind].back().index;
      self_.row(ind) = self_buf.row(self);
      ancestor_.row(ind) = ancestor_buf.row(ancestor);
    }

    self_ = utils::standardise(self_, true, true, false);
    ancestor_ = utils::standardise(ancestor_, true, true, false);

    data_ =
        (self_.transpose() * ancestor_) / static_cast<double>(paths.size() - 1);
  }

 private:
  Component type_;
  std::size_t degree_;
  std::size_t n_ind_;
  std::size_t n_pheno_;
  std::deque<PhenoBuf> history_;
  Eigen::MatrixXd self_;
  Eigen::MatrixXd ancestor_;

  void syncPhenotypes(const State& state);
};

inline void EstimatorAncestorCov::syncPhenotypes(const State& state) {
  if (state.gen == 1) {
    history_.emplace_front(state.pheno(Generation::Parents));
    history_.emplace_front(state.pheno(Generation::Current));
  } else {
    if (history_.size() == degree_ + 1) history_.pop_back();
    history_.emplace_front(state.pheno(Generation::Current));
  }
}

}  // namespace details

inline PopulationEstimator PopulationCousinCov(
    std::size_t degree = 1, Component type = Component::Total) {
  return PopulationEstimator{
      .name = "cousin-" + std::to_string(degree) + "-cov-" + to_string(type),
      .fn = [degree, type](const Params& params) {
        return std::make_unique<details::EstimatorCousinCov>(
            params, degree, type);
      }};
}

inline PopulationEstimator PopulationAncestorCov(
    std::size_t degree = 1, Component type = Component::Total) {
  return PopulationEstimator{
      .name = "ancestor-" + std::to_string(degree) + "-cov-" + to_string(type),
      .fn = [degree, type](const Params& params) {
        return std::make_unique<details::EstimatorAncestorCov>(
            params, degree, type);
      }};
}

}  // namespace amsim
