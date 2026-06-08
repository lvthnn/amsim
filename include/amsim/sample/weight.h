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
