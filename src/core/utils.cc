#include <amsim/core/utils.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>

namespace amsim::utils {
namespace {
void bitmatrix_swap(
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
}  // namespace

void bitmatrix_transpose(std::uint64_t* matrix) {
  std::size_t swap_width = 64;
  auto swap_mask = static_cast<std::uint64_t>(-1);
  while (swap_width != 1) {
    swap_width >>= 1;
    swap_mask = swap_mask ^ (swap_mask >> swap_width);
    bitmatrix_swap(matrix, swap_width, swap_mask);
  }
}

std::vector<std::size_t> order(const Eigen::VectorXd& v) {
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::ranges::stable_sort(
      idx, [&v](std::size_t i0, std::size_t i1) { return v(i0) < v(i1); });
  return idx;
}

Eigen::MatrixXd standardise(const Eigen::MatrixXd& mat) {
  Eigen::RowVectorXd mean = mat.colwise().mean();
  Eigen::RowVectorXd std =
      ((mat.rowwise() - mean).array().square().colwise().sum() / mat.rows())
          .sqrt()
          .max(1e-10);
  return ((mat.rowwise() - mean).array().rowwise() / std.array());
}

}  // namespace amsim::utils
