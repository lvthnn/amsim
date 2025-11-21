#ifndef AMSIMCPP_STATS_H
#define AMSIMCPP_STATS_H

#include <optional>

namespace amsim::stats {

// Sum over the elements of a vector using Kahan summation.
double sum(int N, const double* X, int incX);

// Compute the mean of a vector.
double mean(int N, const double* X, int incX);

// Compute the variance of a vector using Welford's algorithm.
double var(int N, const double* X, int incX, bool population = true);

// Compute the standard deviation of a vector
double std(int N, const double* X, int incX, bool population = true);

// Compute the standard error of a vector
double sem(int N, const double* X, int incX);

// Centre a vector to zero expectation.
void centre(
    int N,
    const double* X,
    int incX,
    double* Y,
    int incY,
    std::optional<double> centre = std::nullopt);

// Scale a vector to unit variance.
void scale(
    int N,
    const double* X,
    int incX,
    double* Y,
    int incY,
    std::optional<double> scale = std::nullopt);

// Standardise a vector.
void standardise(
    int N,
    const double* X,
    int incX,
    double* Y,
    int incY,
    std::optional<double> centre = std::nullopt,
    std::optional<double> scale = std::nullopt);

// Compute the correlation between two vectors
double cor(
    int N,
    const double* X,
    int incX,
    const double* Y,
    int incY,
    std::optional<double> centre_X = std::nullopt,
    std::optional<double> scale_X = std::nullopt,
    std::optional<double> centre_Y = std::nullopt,
    std::optional<double> scale_Y = std::nullopt);

// Compute the cross-correlation matrix of two matrices
void cor(
    int N,
    int P,
    int Q,
    const double* X,
    int ldX,
    const double* Y,
    int ldY,
    double* R,
    int ldR);

// Compute the correlation matrix of a single matrix
void cor(int N, int P, const double* X, int ldX, double* R, int ldR);

double quantile(double q, int N, const double* X, int incX);

}  // namespace amsim::stats

#endif  // AMSIMCPP_STATS_H
