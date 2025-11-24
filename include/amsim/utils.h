#ifndef AMSIMCPP_UTILS_H
#define AMSIMCPP_UTILS_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>

/// Utility functions for simulation operations
namespace amsim::utils {
/// @brief Transpose a 64x64 bit matrix in-place
/// @param matrix Pointer to 64x64 bit matrix
void bitmatrix_transpose(std::uint64_t* matrix);

/// @brief Compute sorting indices for a vector
/// @tparam T Element type
/// @param v Input vector
/// @return Vector of indices that would sort the input
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

/// @brief Assert values are valid probabilities (in [0,1])
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Stride between elements
void assert_probs(std::size_t N, const double* X, std::size_t incX);

/// @brief Assert values are valid correlations (in [-1,1])
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Stride between elements
void assert_cors(std::size_t N, const double* X, std::size_t incX);

/// @brief Assert matrix has unit diagonal
/// @param N Matrix dimension
/// @param X Pointer to matrix
/// @param ldX Leading dimension
void assert_udiag(std::size_t N, const double* X, std::size_t ldX);

/// @brief Assert matrix is positive semi-definite
/// @param N Matrix dimension
/// @param X Pointer to matrix
/// @param ldX Leading dimension
void assert_psd(std::size_t N, const double* X, std::size_t ldX);

/// @brief Assert matrix is a valid correlation matrix
/// @param N Matrix dimension
/// @param X Pointer to matrix
/// @param ldX Leading dimension
void assert_cor(std::size_t N, const double* X, std::size_t ldX);

/// @brief Assert matrix is a valid cross-correlation matrix
/// @param N Matrix dimension
/// @param X Pointer to matrix
/// @param ldX Leading dimension
void assert_cross_cor(std::size_t N, const double* X, std::size_t ldX);

/// @brief Generate random effect sizes
/// @param n_loci Number of loci
void random_effects(std::size_t n_loci);

/// @brief Generate uniform effect sizes
/// @param n_loci Number of loci
void uniform_effects(std::size_t n_loci);

/// @brief Generate random minor allele frequencies
/// @param n_loci Number of loci
void random_mafs(std::size_t n_loci);

/// @brief Generate uniform minor allele frequencies
/// @param n_loci Number of loci
void uniform_mafs(std::size_t n_loci);

}  // namespace amsim::utils

#endif  // AMSIMCPP_UTILS_H
