#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

namespace amsim::utils {

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

inline std::vector<std::string> label_matrix(
    const std::vector<std::string>& row_labels,
    const std::vector<std::string>& col_labels,
    const std::optional<std::string>& row_suffix = std::nullopt,
    const std::optional<std::string>& col_suffix = std::nullopt) {
  std::size_t n_rows = row_labels.size();
  std::size_t n_cols = col_labels.size();
  std::vector<std::string> labels(n_rows * n_cols);

  for (std::size_t row = 0; row < n_rows; ++row)
    for (std::size_t col = 0; col < n_cols; ++col)
      labels[(row * n_cols) + col] =
          (row_labels[row] + row_suffix.value_or("")) +
          "::" + (col_labels[col] + col_suffix.value_or(""));

  return labels;
}

inline std::vector<std::size_t> order(const Eigen::VectorXd& v) {
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::ranges::stable_sort(
      idx, [&v](std::size_t i0, std::size_t i1) { return v(i0) < v(i1); });
  return idx;
}

inline Eigen::MatrixXd standardise(
    const Eigen::MatrixXd& mat, bool population = true) {
  Eigen::RowVectorXd mean = mat.colwise().mean();
  Eigen::RowVectorXd std =
      ((mat.rowwise() - mean).array().square().colwise().sum() /
       (population ? mat.rows() : mat.rows() - 1))
          .sqrt()
          .max(1e-10);
  return ((mat.rowwise() - mean).array().rowwise() / std.array());
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
