#pragma once

#include <amsim/core/distributions.h>

#include <Eigen/Dense>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <filesystem>
#include <fstream>
#include <vector>

namespace amsim {

inline std::vector<std::string> split_string(
    const std::string& s, const std::string& delim = ",") {
  std::vector<std::string> split;
  boost::split(split, s, boost::is_any_of(delim));
  return split;
}

inline Eigen::MatrixXd parse_matrix_value(const std::string& s) {
  std::string delim = (s.find(';') == std::string::npos) ? "\n" : ";";
  std::vector<std::string> rows = split_string(s, delim);
  std::vector<std::vector<double>> matrix(rows.size());

  for (std::size_t row = 0; row < rows.size(); ++row) {
    std::vector<std::string> row_vals = split_string(rows[row]);
    std::ranges::transform(
        row_vals,
        std::back_inserter(matrix[row]),
        [](const std::string& s_val) { return std::stod(s_val); });
  }

  std::size_t n_rows = matrix.size();
  std::size_t n_cols = matrix[0].size();

  if (n_rows == 1) {
    Eigen::VectorXd result(n_cols);
    for (std::size_t col = 0; col < n_cols; ++col) result(col) = matrix[0][col];
    return result;
  }

  Eigen::MatrixXd result(n_rows, n_cols);
  for (std::size_t row = 0; row < n_rows; ++row)
    for (std::size_t col = 0; col < n_cols; ++col)
      result(row, col) = matrix[row][col];

  return result;
}

inline Eigen::MatrixXd parse_matrix_file(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::string contents(
      (std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
  return parse_matrix_value(contents);
}

inline std::vector<std::size_t> parse_indices(const std::string& s) {
  std::vector<std::string> tokens = split_string(s);
  std::vector<std::size_t> indices(tokens.size());
  std::ranges::transform(tokens, indices.begin(), [](const std::string& t) {
    return static_cast<std::size_t>(std::stoull(t));
  });
  return indices;
}

inline std::vector<std::size_t> parse_indices_file(
    const std::filesystem::path& path) {
  std::ifstream file(path);
  std::string contents(
      (std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
  return parse_indices(contents);
}

template <typename T>
struct File {
  std::filesystem::path path;
  T load() const {
    if constexpr (std::is_same_v<T, std::vector<std::size_t>>)
      return parse_indices_file(path);
    else if constexpr (std::is_same_v<T, Eigen::MatrixXd>)
      return parse_matrix_file(path);
    else
      static_assert(false, "Unsupported File type");
  }
};


inline amsim::Distribution parse_distribution(
    const std::string& s, bool is_probability = false) {
  int paren_begin = s.find('(');
  int paren_end = s.find(')');
  if (paren_begin == std::string::npos || paren_end == std::string::npos ||
      paren_begin >= paren_end)
    throw std::runtime_error("Invalid distribution form");

  std::string dist_name = s.substr(0, paren_begin);
  std::string params_str =
      s.substr(paren_begin + 1, paren_end - paren_begin - 1);
  std::vector<std::string> params_str_vec = split_string(params_str);
  std::vector<double> params(params_str_vec.size());

  if (!params_str.empty())
    std::ranges::transform(
        params_str_vec, params.begin(), [](const std::string& s) {
          return std::stod(s);
        });

  boost::to_lower(dist_name);

  return amsim::str_to_distribution(dist_name, params, is_probability);
}

}  // namespace amsim
