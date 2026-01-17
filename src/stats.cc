#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

#include <amsim/stats.h>

namespace amsim::stats {

double sum(int N, const double* X, int incX) {
  double y;
  double t;
  double s = 0.0;
  double c = 0.0;
  for (int i = 0; i < N; ++i, X += incX) {
    y = *X - c;
    t = s + y;
    c = (t - s) - y;
    s = t;
  }
  return s;
}

double mean(int N, const double* X, int incX) {
  return sum(N, X, incX) / static_cast<double>(N);
}

double var(int N, const double* X, int incX, bool population) {
  double mean = 0.0;
  double m2 = 0.0;
  double x;
  double delta;
  double delta2;

  const double* xi = X;

  for (int i = 0; i < N; ++i, xi += incX) {
    x = *xi;
    delta = x - mean;
    mean += delta / (i + 1);
    delta2 = x - mean;
    m2 += delta * delta2;
  }

  double denom =
      population ? static_cast<double>(N) : static_cast<double>(N - 1);
  return m2 / denom;
}

double std(int N, const double* X, int incX, bool population) {
  return std::sqrt(var(N, X, incX, population));
}

double sem(int N, const double* X, int incX) {
  return std(N, X, incX, false) / std::sqrt(static_cast<double>(N));
}

void centre(
    int N,
    const double* X,
    int incX,
    double* Y,
    int incY,
    std::optional<double> centre) {
  if (!centre) centre = mean(N, X, incX);
  for (int i = 0; i < N; ++i) Y[i * incY] = X[i * incX] - *centre;
}

void scale(
    int N,
    const double* X,
    int incX,
    double* Y,
    int incY,
    std::optional<double> scale) {
  if (!scale) scale = std::sqrt(var(N, X, incX));
  const double* xi = X;
  double* eta = Y;
  for (int i = 0; i < N; ++i, xi += incX, eta += incY) *eta = *xi / *scale;
}

void standardise(
    int N,
    const double* X,
    int incX,
    double* Y,
    int incY,
    std::optional<double> centre_val,
    std::optional<double> scale_val) {
  double c = centre_val.value_or(mean(N, X, incX));
  centre(N, X, incX, Y, incY, c);

  double s = scale_val.value_or(std::sqrt(var(N, Y, incY, true)));
  if (s > 0.0) {
    cblas_dscal(N, 1.0 / s, Y, incY);
  }
}

double cor(
    int N,
    const double* X,
    int incX,
    const double* Y,
    int incY,
    std::optional<double> centre_X,
    std::optional<double> scale_X,
    std::optional<double> centre_Y,
    std::optional<double> scale_Y) {
  if (!centre_X) centre_X = mean(N, X, incX);
  if (!centre_Y) centre_Y = mean(N, Y, incY);
  if (!scale_X) scale_X = std::sqrt(var(N, X, incX));
  if (!scale_Y) scale_Y = std::sqrt(var(N, Y, incY));
  double cor = 0.0;
  double denom = 1.0 / (static_cast<double>(N) * *scale_X * *scale_Y);
  const double* xi = X;
  const double* eta = Y;

  for (int i = 0; i < N; ++i, xi += incX, eta += incY)
    cor += denom * (*xi - *centre_X) * (*eta - *centre_Y);

  return cor;
}

void cor(
    int N,
    int P,
    int Q,
    const double* X,
    int ldX,
    const double* Y,
    int ldY,
    double* R,
    int ldR) {
  auto dn = static_cast<double>(N);

  // compute column sums of X and Y (sx and sy)
  std::vector<double> ones(N, 1.0);
  std::vector<double> sx(P);
  std::vector<double> sy(Q);

  cblas_dgemv(
      CblasColMajor,
      CblasTrans,
      N,
      P,
      1.0,
      X,
      ldX,
      ones.data(),
      1,
      0.0,
      sx.data(),
      1);

  cblas_dgemv(
      CblasColMajor,
      CblasTrans,
      N,
      Q,
      1.0,
      Y,
      ldY,
      ones.data(),
      1,
      0.0,
      sy.data(),
      1);

  // compute column sum of squares of X and Y (sxx and syy)
  std::vector<double> sxx(P);
  std::vector<double> syy(Q);
  double nrm_i;
  double nrm_j;

  for (int i = 0; i < P; ++i) {
    nrm_i = cblas_dnrm2(N, X + (i * ldX), 1);
    sxx[i] = nrm_i * nrm_i;
  }

  for (int j = 0; j < Q; ++j) {
    nrm_j = cblas_dnrm2(N, Y + (j * ldY), 1);
    syy[j] = nrm_j * nrm_j;
  }

  // compute X^TY and store in buffer R
  cblas_dgemm(
      CblasColMajor,
      CblasTrans,
      CblasNoTrans,
      P,
      Q,
      N,
      1.0,
      X,
      ldX,
      Y,
      ldY,
      0.0,
      R,
      ldR);

  // subtract rank-one 1/n * sx * sy^T
  cblas_dger(
      CblasColMajor, P, Q, -1.0 / dn, sx.data(), 1, sy.data(), 1, R, ldR);

  // compute standard deviations
  std::vector<double> sdx(P);
  std::vector<double> sdy(Q);

  for (int j = 0; j < P; ++j) {
    double vx = sxx[j] - ((sx[j] * sx[j]) / dn);
    sdx[j] = std::sqrt(std::max(0.0, vx));
  }

  for (int k = 0; k < Q; ++k) {
    double vy = syy[k] - ((sy[k] * sy[k]) / dn);
    sdy[k] = std::sqrt(std::max(0.0, vy));
  }

  // divide entries by standard deviations
  for (int i = 0; i < P; ++i) {
    if (sdx[i] == 0.0)
      for (int j = 0; j < Q; ++j)
        R[i + (j * ldR)] = std::numeric_limits<double>::quiet_NaN();
    else
      cblas_dscal(Q, 1.0 / sdx[i], R + i, ldR);
  }
  for (int j = 0; j < Q; ++j) {
    if (sdy[j] == 0.0)
      for (int i = 0; i < P; ++i)
        R[i + (j * ldR)] = std::numeric_limits<double>::quiet_NaN();
    else
      cblas_dscal(P, 1.0 / sdy[j], R + (j * ldR), 1);
  }
}

void cor(int N, int P, const double* X, int ldX, double* R, int ldR) {
  cor(N, P, P, X, ldX, X, ldX, R, ldR);
}

double quantile(const double q, int N, const double* X, int incX) {
  if (q < 0.0) throw std::invalid_argument("can't specify quantile below zero");
  if (q > 1.0) throw std::invalid_argument("can't specify quantile above one");
  if (N <= 0) return std::numeric_limits<double>::quiet_NaN();
  if (N == 1) return X[0];

  // copy into temporary buffer
  std::vector<double> tmp(N);
  const double* xi = X;
  for (int i = 0; i < N; ++i, xi += incX) tmp[i] = *xi;

  // R type 7 quantile interpolation strategy
  const double h = 1.0 + ((N - 1) * q);
  int lo = static_cast<int>(std::floor(h)) - 1;
  int hi = static_cast<int>(std::ceil(h)) - 1;
  lo = std::max(0, lo);
  hi = std::max(0, hi);
  if (lo >= N) lo = N - 1;
  if (hi >= N) hi = N - 1;

  const double gamma = std::max(0.0, h - std::floor(h));

  std::nth_element(tmp.begin(), tmp.begin() + lo, tmp.end());
  const double x_lo = tmp[lo];

  if (gamma == 0.0 || hi == lo) return x_lo;

  std::nth_element(tmp.begin(), tmp.begin() + hi, tmp.end());
  const double x_hi = tmp[hi];

  return ((1.0 - gamma) * x_lo) + (gamma * x_hi);
}

}  // namespace amsim::stats
