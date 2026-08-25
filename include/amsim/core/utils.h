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
#include <sys/wait.h>

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <semaphore>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace amsim::utils {

inline std::vector<std::string> split_string(
    const std::string& s, const std::string& delim = ",") {
  std::vector<std::string> split;
  boost::split(split, s, boost::is_any_of(delim));
  std::erase_if(split, [](const std::string& s) { return s.empty(); });
  return split;
}

inline void bitmatrix_swap(
    std::uint64_t matrix[], std::size_t width, std::uint64_t mask) {
  std::size_t inner;
  std::size_t outer;
  for (outer = 0; outer < 64 / (width * 2); ++outer) {
    for (inner = 0; inner < width; ++inner) {
      std::uint64_t* x = &matrix[(inner) + (outer * width * 2)];
      std::uint64_t* y = &matrix[(inner + width) + (outer * width * 2)];
      *x = ((*y << width) & mask) ^ *x;
      *y = ((*x & mask) >> width) ^ *y;
      *x = ((*y << width) & mask) ^ *x;
    }
  }
}

inline void bitmatrix_transpose(std::uint64_t* matrix) {
  std::size_t swap_width = 64;
  auto swap_mask = static_cast<std::uint64_t>(-1);
  while (swap_width != 1) {
    swap_width >>= 1;
    swap_mask = swap_mask ^ (swap_mask >> swap_width);
    bitmatrix_swap(matrix, swap_width, swap_mask);
  }
}

inline Eigen::MatrixXd random_orthogonal(std::size_t dim) {
  Eigen::MatrixXd random = Eigen::MatrixXd::Random(dim, dim);
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(random);
  return qr.householderQ();
}

inline Eigen::MatrixXd matrix_from_singular_values(
    const std::string& s, bool symmetric = false) {
  std::vector<std::string> vs = split_string(s);
  std::size_t n_pheno = vs.size();
  Eigen::MatrixXd u_mat;
  Eigen::MatrixXd v_mat;

  Eigen::VectorXd singular_values(n_pheno);
  std::ranges::transform(
      vs, singular_values.begin(), [](const std::string& s_val) {
        return std::stod(s_val);
      });

  Eigen::MatrixXd s_mat = singular_values.asDiagonal();

  u_mat = random_orthogonal(n_pheno);
  if (!symmetric) v_mat = random_orthogonal(n_pheno);

  return (symmetric) ? u_mat * s_mat * u_mat.transpose()
                     : u_mat * s_mat * v_mat.transpose();
}

inline Eigen::MatrixXd standardise(
    const Eigen::MatrixXd& mat,
    bool population = true,
    bool centre = true,
    bool scale = true) {
  Eigen::RowVectorXd mean = centre ? Eigen::RowVectorXd(mat.colwise().mean())
                                   : Eigen::RowVectorXd::Zero(mat.cols());
  if (!scale) return mat.rowwise() - mean;
  Eigen::RowVectorXd std =
      ((mat.rowwise() - mean).array().square().colwise().sum() /
       (population ? mat.rows() : mat.rows() - 1))
          .sqrt()
          .max(1e-10);
  return ((mat.rowwise() - mean).array().rowwise() / std.array());
}

inline std::vector<std::string> vector_prefix(
    const std::vector<std::string>& labels,
    const std::optional<std::string>& prefix = std::nullopt) {
  std::vector<std::string> result(labels.size());
  for (std::size_t i = 0; i < labels.size(); ++i)
    result[i] = prefix.value_or("") + labels[i];
  return result;
}

inline std::vector<std::string> vector_suffix(
    const std::vector<std::string>& labels,
    const std::optional<std::string>& suffix = std::nullopt) {
  std::vector<std::string> result(labels.size());
  for (std::size_t i = 0; i < labels.size(); ++i)
    result[i] = labels[i] + suffix.value_or("");
  return result;
}

inline std::vector<std::size_t> order(const Eigen::VectorXd& v) {
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::ranges::stable_sort(
      idx, [&v](std::size_t i0, std::size_t i1) { return v(i0) < v(i1); });
  return idx;
}

inline std::counting_semaphore<64> process_semaphore{
    static_cast<std::ptrdiff_t>(
        std::thread::hardware_concurrency() > 0
            ? std::thread::hardware_concurrency()
            : 4)};

inline void system_throttled(const std::string& cmd) {
  process_semaphore.acquire();
  FILE* pipe = popen((cmd + " 2>&1").c_str(), "r");
  std::string output;
  if (pipe) {
    char buffer[512];
    while (fgets(buffer, sizeof(buffer), pipe)) output += buffer;
  }
  int rc = pclose(pipe);
  process_semaphore.release();
  if (WIFEXITED(rc) && WEXITSTATUS(rc) != 0) {
    throw std::runtime_error(output);
  }
}

inline void check_plink2() {
  int rc = std::system("command -v plink2 >/dev/null 2>&1");
  if (!WIFEXITED(rc) || WEXITSTATUS(rc) != 0)
    throw std::runtime_error("Executable 'plink2' not found in PATH");
}

inline void check_gcta64() {
  int rc = std::system("command -v gcta64 >/dev/null 2>&1");
  if (!WIFEXITED(rc) || WEXITSTATUS(rc) != 0)
    throw std::runtime_error("Executable 'gcta64' not found in PATH");
}

}  // namespace amsim::utils
