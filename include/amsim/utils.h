#ifndef AMSIMCPP_UTILS_H
#define AMSIMCPP_UTILS_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>

namespace amsim::utils {
void bitmatrix_transpose(std::uint64_t* matrix);

template <typename T>
std::vector<std::size_t> order(const std::vector<T>& v) {
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(
      idx.begin(), idx.end(), [&v](std::size_t i0, std::size_t i1) {
        return v[i0] < v[i1];
      });
  return idx;
}

void assert_probs(std::size_t N, const double* X, std::size_t incX);
void assert_cors(std::size_t N, const double* X, std::size_t incX);
void assert_udiag(std::size_t N, const double* X, std::size_t ldX);
void assert_psd(std::size_t N, const double* X, std::size_t ldX);
void assert_cor(std::size_t N, const double* X, std::size_t ldX);
void assert_cross_cor(std::size_t N, const double* X, std::size_t ldX);

void random_effects(std::size_t n_loci);
void uniform_effects(std::size_t n_loci);
void random_mafs(std::size_t n_loci);
void uniform_mafs(std::size_t n_loci);

}  // namespace amsim::utils

#endif  // AMSIMCPP_UTILS_H
