#ifndef AMSIMCPP_STATS_H
#define AMSIMCPP_STATS_H

#pragma once

#include <optional>

namespace amsim::stats {

// Sum over the elements of a vector using Kahan summation.
double sum(const int N, const double* X, const int incX);

// Compute the mean of a vector.
double mean(const int N, const double* X, const int incX);

// Compute the variance of a vector using Welford's algorithm.
double var(
    const int N, const double* X, const int incX, bool population = true);

// Compute the standard deviation of a vector
double std(
    const int N, const double* X, const int incX, bool population = true);

// Compute the standard error of a vector
double sem(const int N, const double* X, const int incX);

// Centre a vector to zero expectation.
void centre(
    const int N,
    const double* X,
    const int incX,
    double* Y,
    const int incY,
    std::optional<double> centre = std::nullopt);

// Scale a vector to unit variance.
void scale(
    const int N,
    const double* X,
    const int incX,
    double* Y,
    const int incY,
    std::optional<double> scale = std::nullopt);

// Standardise a vector.
void standardise(
    const int N,
    const double* X,
    const int incX,
    double* Y,
    const int incY,
    std::optional<double> centre = std::nullopt,
    std::optional<double> scale = std::nullopt);

// Compute the correlation between two vectors
double cor(
    const int N,
    const double* X,
    const int incX,
    const double* Y,
    const int incY,
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
    const int ldX,
    const double* Y,
    const int ldY,
    double* R,
    const int ldR);

// Compute the correlation matrix of a single matrix
void cor(
    int N, int P, const double* X, const int ldX, double* R, const int ldR);

double quantile(double q, int N, const double* X, const int ldX);

}  // namespace amsim::stats

#endif  // AMSIMCPP_STATS_H
