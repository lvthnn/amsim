#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

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

void assert_probs(std::size_t N, const double* X, std::size_t incX) {
  const double* xi = X;

  for (std::size_t i = 0; i < N; ++i, xi += incX)
    if (*xi < 0.0 || *xi > 1.0)
      throw std::invalid_argument(
          "invalid element " + std::to_string(i) + "; must be in [0, 1]");
}

void assert_cors(std::size_t N, const double* X, std::size_t incX) {
  const double* xi = X;

  for (std::size_t i = 0; i < N; ++i, xi += incX)
    if (*xi < -1.0 || *xi > 1.0)
      throw std::invalid_argument(
          "invalid element " + std::to_string(i) + "; must be in [-1, 1]");
}

void assert_udiag(std::size_t N, const double* X, std::size_t ldX) {
  const double* xi = X;

  for (std::size_t i = 0; i < N; ++i, xi += ldX) {
    const double* el_diag = xi + i;
    if (std::abs(*el_diag - 1.0) > 1e-12)
      throw std::invalid_argument(
          "diagonal element " + std::to_string(i) + " does not equal one");
  }
}

void assert_psd(std::size_t N, const double* X, std::size_t ldX) {
  std::vector<double> x_copy(N * N);
  std::vector<double> eigvals(N);
  std::copy_n(X, N * N, x_copy.begin());

  const char clpk_job = 'N';
  const char clpk_uplo = 'L';
  const int clpk_ld = static_cast<int>(ldX);
  int clpk_info;

  double clpk_wkopt;
  int clpk_lwork = -1;

#ifdef __APPLE__
  dsyev_(
      &clpk_job,
      &clpk_uplo,
      &clpk_ld,
      x_copy.data(),
      &clpk_ld,
      eigvals.data(),
      &clpk_wkopt,
      &clpk_lwork,
      &clpk_info);
#else
  LAPACK_dsyev(
      &clpk_job,
      &clpk_uplo,
      &clpk_ld,
      X_copy.data(),
      &clpk_ld,
      eigvals.data(),
      &clpk_work_query,
      &clpk_lwork,
      &clpk_info);
#endif

  clpk_lwork = static_cast<int>(clpk_wkopt);
  std::vector<double> clpk_work(clpk_lwork);

#ifdef __APPLE__
  dsyev_(
      &clpk_job,
      &clpk_uplo,
      &clpk_ld,
      x_copy.data(),
      &clpk_ld,
      eigvals.data(),
      clpk_work.data(),
      &clpk_lwork,
      &clpk_info);
#else
  LAPACK_dsyev(
      &clpk_job,
      &clpk_uplo,
      &clpk_ld,
      X_copy.data(),
      &clpk_ld,
      eigvals.data(),
      clpk_work.data(),
      &clpk_lwork,
      &clpk_info);
#endif

  if (clpk_info != 0)
    throw std::runtime_error(
        "LAPACK dsyev failed with INFO = " + std::to_string(clpk_info));

  for (const double& eigval : eigvals)
    if (eigval < -1e-12)
      throw std::invalid_argument("matrix not positive semi-definite");
}

void assert_cor(std::size_t N, const double* X, std::size_t ldX) {
  assert_cors(N * N, X, 1);
  assert_udiag(N, X, ldX);
  assert_psd(N, X, ldX);
}

void assert_cross_cor(std::size_t N, const double* X, std::size_t ldX) {
  assert_cors(N * N, X, 1);

  std::vector<double> x_copy(N * N);
  std::vector<double> svals(N);
  std::copy_n(X, N * N, x_copy.begin());

  const char clpk_job = 'N';
  const int clpk_ld = static_cast<int>(ldX);
  int clpk_lwork = -1;
  int clpk_info;
  double clpk_wkopt;

#ifdef __APPLE__
  dgesvd_(
      &clpk_job,
      &clpk_job,
      &clpk_ld,
      &clpk_ld,
      x_copy.data(),
      &clpk_ld,
      svals.data(),
      nullptr,
      &clpk_ld,
      nullptr,
      &clpk_ld,
      &clpk_wkopt,
      &clpk_lwork,
      &clpk_info);
#else
  LAPACK_dgesvd(
      &clpk_job,
      &clpk_job,
      &clpk_ld,
      &clpk_ld,
      X_copy.data(),
      &clpk_ld,
      svals.data(),
      nullptr,
      &clpk_ld,
      nullptr,
      &clpk_ld,
      &clpk_wkopt,
      &clpk_lwork,
      &clpk_info);
#endif

  if (clpk_info != 0) {
    throw std::runtime_error(
        "lapack dgesvd failed with info = " + std::to_string(clpk_info));
  }

  clpk_lwork = static_cast<int>(clpk_wkopt);
  std::vector<double> clpk_work(clpk_lwork);

#ifdef __APPLE__
  dgesvd_(
      &clpk_job,
      &clpk_job,
      &clpk_ld,
      &clpk_ld,
      x_copy.data(),
      &clpk_ld,
      svals.data(),
      nullptr,
      &clpk_ld,
      nullptr,
      &clpk_ld,
      clpk_work.data(),
      &clpk_lwork,
      &clpk_info);
#else
  LAPACK_dgesvd(
      &clpk_job,
      &clpk_job,
      &clpk_ld,
      &clpk_ld,
      X_copy.data(),
      &clpk_ld,
      svals.data(),
      nullptr,
      &clpk_ld,
      nullptr,
      &clpk_ld,
      clpk_work.data(),
      &clpk_lwork,
      &clpk_info);
#endif

  if (clpk_info != 0)
    throw std::runtime_error(
        "LAPACK dgesvd failed with INFO = " + std::to_string(clpk_info));

  for (const double& sval : svals) {
    if (sval > 1.0)
      throw std::invalid_argument(
          "infeasible between-mate correlation regime; singular value exceeds "
          "one (" +
          std::to_string(sval) + ")");
  }
}

// function to generate evenly-spaced chromosome separation sites for
// recombination

// function to generate UNLINKED recombination map with constant MAFs

// function to generate UNLINKED recombination map with random MAFs

// function to generate LINKED recombination map with otherwise constant MAFs
// -> fill in the indicated chromosome separation sites

// function to generate LINKED recombination map with random MAFs
// -> use beta? We don't have that one...

// function to generate phenotype effect vector with uniform effect sizes
// -> constant vector of heritability divided by sqrt of no. loci

// function to generate phenotype effect vector with random effect sizes
// -> use polar normal RNG, normalise, and scale

}  // namespace amsim::utils
