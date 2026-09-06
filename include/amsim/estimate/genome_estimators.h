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
#include <amsim/estimate/estimator.h>
#include <amsim/estimate/strategy.h>

#include <Eigen/Dense>

namespace amsim {

namespace details {

class EstimatorGenotypeMean : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorGenotypeMean(const Params& params)
      : PopulationEstimatorStrategy("geno_mean", {}, {}, params.geno.n_loc) {};

  void compute(const State& state) override { data_ = state.geno().v_lmean(); }
};

class EstimatorGenotypeVar : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorGenotypeVar(const Params& params)
      : PopulationEstimatorStrategy("geno_var", {}, {}, params.geno.n_loc) {}

  void compute(const State& state) override { data_ = state.geno().v_lvar(); }
};

class EstimatorGenotypeMAF : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorGenotypeMAF(const Params& params)
      : PopulationEstimatorStrategy("geno_maf", {}, {}, params.geno.n_loc) {}

  void compute(const State& state) override { data_ = state.geno().v_lmaf(); }
};

class EstimatorGenotypeCov : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorGenotypeCov(const Params& params)
      : PopulationEstimatorStrategy(
            "geno_cov", {}, {}, params.geno.n_loc, params.geno.n_loc) {
    n_loc_ = params.geno.n_loc;
  }

  void compute(const State& state) override {
    const GenoBuf& geno = state.geno();
    std::size_t n_words = geno.n_words();
    std::size_t n_ind = geno.n_ind();

    for (std::size_t loc1 = 0; loc1 < n_loc_; ++loc1) {
      // haplotypes of the first locus
      const uint64_t* h01 = geno.h0().rowptr(loc1);
      const uint64_t* h11 = geno.h1().rowptr(loc1);
      for (std::size_t loc2 = loc1; loc2 < n_loc_; ++loc2) {
        // haplotypes of the second locus
        const uint64_t* h02 = geno.h0().rowptr(loc2);
        const uint64_t* h12 = geno.h1().rowptr(loc2);
        std::size_t acc = 0;

        for (std::size_t word = 0; word < n_words; ++word)
          if (word == n_words - 1 && (n_ind % 64)) {
            std::uint64_t mask = (1ULL << (n_ind % 64)) - 1ULL;
            acc += __builtin_popcountll(h01[word] & h02[word] & mask) +
                   __builtin_popcountll(h01[word] & h12[word] & mask) +
                   __builtin_popcountll(h11[word] & h02[word] & mask) +
                   __builtin_popcountll(h11[word] & h12[word] & mask);
          } else {
            acc += __builtin_popcountll(h01[word] & h02[word]) +
                   __builtin_popcountll(h01[word] & h12[word]) +
                   __builtin_popcountll(h11[word] & h02[word]) +
                   __builtin_popcountll(h11[word] & h12[word]);
          }

        data_(loc1, loc2) =
            ((1.0 / n_ind) * acc) - (geno.v_lmean(loc1) * geno.v_lmean(loc2));
      }
    }

    data_.triangularView<Eigen::Lower>() =
        data_.transpose().triangularView<Eigen::Lower>();
  }

 private:
  std::size_t n_loc_;
};

class EstimatorGenotypeCor : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorGenotypeCor(const Params& params)
      : PopulationEstimatorStrategy(
            "geno_cor", {}, {}, params.geno.n_loc, params.geno.n_loc) {
    n_loc_ = params.geno.n_loc;
  }

  void compute(const State& state) override {
    const GenoBuf& geno = state.geno();
    std::size_t n_words = geno.n_words();
    std::size_t n_ind = geno.n_ind();

    for (std::size_t loc1 = 0; loc1 < n_loc_; ++loc1) {
      const uint64_t* h01 = geno.h0().rowptr(loc1);
      const uint64_t* h11 = geno.h1().rowptr(loc1);
      for (std::size_t loc2 = loc1; loc2 < n_loc_; ++loc2) {
        const uint64_t* h02 = geno.h0().rowptr(loc2);
        const uint64_t* h12 = geno.h1().rowptr(loc2);
        std::size_t acc = 0;

        for (std::size_t word = 0; word < n_words; ++word) {
          if (word == n_words - 1 && (n_ind % 64)) {
            std::uint64_t mask = (1ULL << (n_ind % 64)) - 1ULL;
            acc += __builtin_popcountll(h01[word] & h02[word] & mask) +
                   __builtin_popcountll(h01[word] & h12[word] & mask) +
                   __builtin_popcountll(h11[word] & h02[word] & mask) +
                   __builtin_popcountll(h11[word] & h12[word] & mask);
          } else {
            acc += __builtin_popcountll(h01[word] & h02[word]) +
                   __builtin_popcountll(h01[word] & h12[word]) +
                   __builtin_popcountll(h11[word] & h02[word]) +
                   __builtin_popcountll(h11[word] & h12[word]);
          }
        }

        double cov =
            ((1.0 / n_ind) * acc) - (geno.v_lmean(loc1) * geno.v_lmean(loc2));
        double denom = std::sqrt(geno.v_lvar(loc1) * geno.v_lvar(loc2));
        data_(loc1, loc2) = (denom > 0.0) ? cov / denom : 0.0;
      }
    }

    data_.triangularView<Eigen::Lower>() =
        data_.transpose().triangularView<Eigen::Lower>();
  }

 private:
  std::size_t n_loc_;
};

}  // namespace details

inline PopulationEstimator PopulationGenotypeCov() {
  return PopulationEstimator{
      .name = "genotype-cov", .fn = [](const Params& params) {
        return std::make_unique<details::EstimatorGenotypeCov>(params);
      }};
}

inline PopulationEstimator PopulationGenotypeMean() {
  return PopulationEstimator{
      .name = "genotype-mean", .fn = [](const Params& params) {
        return std::make_unique<details::EstimatorGenotypeMean>(params);
      }};
}

inline PopulationEstimator PopulationGenotypeVar() {
  return PopulationEstimator{
      .name = "genotype-var", .fn = [](const Params& params) {
        return std::make_unique<details::EstimatorGenotypeVar>(params);
      }};
}

inline PopulationEstimator PopulationGenotypeMAF() {
  return PopulationEstimator{
      .name = "genotype-maf", .fn = [](const Params& params) {
        return std::make_unique<details::EstimatorGenotypeMAF>(params);
      }};
}

inline PopulationEstimator PopulationGenotypeCor() {
  return PopulationEstimator{
      .name = "genotype-cor", .fn = [](const Params& params) {
        return std::make_unique<details::EstimatorGenotypeCor>(params);
      }};
}

}  // namespace amsim
