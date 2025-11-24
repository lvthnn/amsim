#ifndef AMSIMCPP_STATS_H
#define AMSIMCPP_STATS_H

#include <optional>

/// Statistical functions for vectors and matrices
namespace amsim::stats {

/// @brief Sum vector elements using Kahan summation
///
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Stride between elements
///
/// @return Sum of elements
double sum(int N, const double* X, int incX);

/// @brief Compute vector mean
///
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Stride between elements
///
/// @return Mean value
double mean(int N, const double* X, int incX);

/// @brief Compute vector variance using Welford's algorithm
///
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Stride between elements
/// @param population Whether to use population variance (default: true)
///
/// @return Variance
double var(int N, const double* X, int incX, bool population = true);

/// @brief Compute vector standard deviation
///
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Stride between elements
/// @param population Whether to use population variance (default: true)
///
/// @return Standard deviation
double std(int N, const double* X, int incX, bool population = true);

/// @brief Compute vector standard error of the mean
///
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Stride between elements
///
/// @return Standard error
double sem(int N, const double* X, int incX);

/// @brief Centre vector to zero mean
///
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Input stride
/// @param Y Pointer to output vector
/// @param incY Output stride
/// @param centre Optional custom centre value
void centre(
    int N,
    const double* X,
    int incX,
    double* Y,
    int incY,
    std::optional<double> centre = std::nullopt);

/// @brief Scale vector to unit variance
///
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Input stride
/// @param Y Pointer to output vector
/// @param incY Output stride
/// @param scale Optional custom scale value
void scale(
    int N,
    const double* X,
    int incX,
    double* Y,
    int incY,
    std::optional<double> scale = std::nullopt);

/// @brief Standardise vector to zero mean and unit variance
///
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Input stride
/// @param Y Pointer to output vector
/// @param incY Output stride
/// @param centre Optional custom centre value
///
/// @param scale Optional custom scale value
void standardise(
    int N,
    const double* X,
    int incX,
    double* Y,
    int incY,
    std::optional<double> centre = std::nullopt,
    std::optional<double> scale = std::nullopt);

/// @brief Compute correlation between two vectors
///
/// @param N Number of elements
/// @param X Pointer to first vector
/// @param incX First vector stride
/// @param Y Pointer to second vector
/// @param incY Second vector stride
/// @param centre_X Optional custom centre for X
/// @param scale_X Optional custom scale for X
/// @param centre_Y Optional custom centre for Y
/// @param scale_Y Optional custom scale for Y
///
/// @return Correlation coefficient
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

/// @brief Compute cross-correlation matrix between two matrices
///
/// @param N Number of rows
/// @param P Number of columns in X
/// @param Q Number of columns in Y
/// @param X Pointer to first matrix
/// @param ldX Leading dimension of X
/// @param Y Pointer to second matrix
/// @param ldY Leading dimension of Y
/// @param R Pointer to output correlation matrix
/// @param ldR Leading dimension of R
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

/// @brief Compute correlation matrix of a single matrix
///
/// @param N Number of rows
/// @param P Number of columns
/// @param X Pointer to input matrix
/// @param ldX Leading dimension of X
/// @param R Pointer to output correlation matrix
/// @param ldR Leading dimension of R
void cor(int N, int P, const double* X, int ldX, double* R, int ldR);

/// @brief Compute quantile of a vector
///
/// @param q Quantile to compute (0 to 1)
/// @param N Number of elements
/// @param X Pointer to input vector
/// @param incX Stride between elements
///
/// @return Quantile value
double quantile(double q, int N, const double* X, int incX);

}  // namespace amsim::stats

#endif  // AMSIMCPP_STATS_H
