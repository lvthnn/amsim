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

#include <amsim/core/state.h>
#include <amsim/io/writer.h>

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace amsim {

class PopulationEstimatorStrategy {
 public:
  PopulationEstimatorStrategy(
      std::string name,
      std::vector<std::string> row_labels,
      std::vector<std::string> col_labels,
      std::size_t n_rows,
      std::size_t n_cols = 1)
      : name_(std::move(name)),
        row_labels_(std::move(row_labels)),
        col_labels_(std::move(col_labels)),
        n_rows_(n_rows),
        n_cols_(n_cols),
        data_(n_rows_, n_cols_) {
    Writer::create(data_, name_, row_labels_, col_labels_);
  }

  virtual ~PopulationEstimatorStrategy() = default;
  virtual void compute(const State& state) = 0;

  std::string name() const { return name_; }
  std::vector<std::string> row_labels() const { return row_labels_; }
  std::vector<std::string> col_labels() const { return col_labels_; }
  std::size_t n_rows() const { return n_rows_; }
  std::size_t n_cols() const { return n_cols_; }

  void operator()(const State& state) {
    compute(state);
    Writer::write(data_, name(), state.rep, state.gen);
  }

 protected:
  std::string name_;
  std::vector<std::string> row_labels_;
  std::vector<std::string> col_labels_;
  std::size_t n_rows_;
  std::size_t n_cols_;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> data_;
};

}  // namespace amsim
