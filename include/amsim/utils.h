#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace amsim::utils {

void bitmatrix_transpose(std::uint64_t* matrix);

std::vector<std::size_t> order(const Eigen::VectorXd& v);

Eigen::MatrixXd standardise(const Eigen::MatrixXd& mat);

Eigen::MatrixXd cor(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y);

}  // namespace amsim::utils
