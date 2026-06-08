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

#include <functional>

#include <Eigen/Dense>

namespace amsim {

using WeightFunction =
    std::function<void(const Eigen::MatrixXd&, Eigen::VectorXd&)>;

inline WeightFunction Uniform() {
  return [](const Eigen::MatrixXd& /*agg*/, Eigen::VectorXd& res) {
    res.setConstant(1.0);
  };
}

inline WeightFunction Logistic(const Eigen::VectorXd& effects) {
  return [effects](const Eigen::MatrixXd& agg, Eigen::VectorXd& res) {
    res = 1.0 / (1.0 + (-agg * effects).array().exp());
  };
}

};  // namespace amsim
