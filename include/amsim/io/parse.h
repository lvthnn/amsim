// This file is part of amsim, copyright (C) 2025-2026 Kári Hlynsson.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <amsim/core/distributions.h>
#include <amsim/sample/weight.h>
#include <amsim/sample/proband.h>
#include <amsim/core/utils.h>

#include <Eigen/Dense>
#include <algorithm>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <filesystem>
#include <fstream>
#include <vector>

namespace amsim {

inline Eigen::MatrixXd parse_matrix_value(const std::string& s) {
  std::string norm = s;
  std::ranges::replace(norm, ';', '\n');
  std::ranges::replace(norm, ',', ' ');

  std::vector<std::string> rows = utils::split_string(norm, "\n");
  std::vector<std::vector<double>> matrix(rows.size());

  for (std::size_t row = 0; row < rows.size(); ++row) {
    std::vector<std::string> row_vals = utils::split_string(rows[row], " ");
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

inline Eigen::MatrixXd parse_pheno_file(const std::filesystem::path& path) {
  std::ifstream file(path);
  if (!file.is_open())
    throw std::runtime_error(
        "Could not open file " + path.string() + " for parsing");

  std::string line;
  std::getline(file, line);  // skip header

  std::vector<std::vector<double>> matrix;
  while (std::getline(file, line)) {
    std::vector<std::string> tokens = utils::split_string(line, "\t");
    std::vector<double> row;
    for (std::size_t i = 2; i < tokens.size(); ++i)
      row.push_back(std::stod(tokens[i]));
    matrix.push_back(std::move(row));
  }

  std::size_t n_rows = matrix.size();
  std::size_t n_cols = matrix[0].size();
  Eigen::MatrixXd result(n_rows, n_cols);
  for (std::size_t r = 0; r < n_rows; ++r)
    for (std::size_t c = 0; c < n_cols; ++c)
      result(r, c) = matrix[r][c];
  return result;
}

inline Eigen::MatrixXd parse_matrix_file(const std::filesystem::path& path) {
  std::ifstream file(path);
  if (!file.is_open())
    throw std::runtime_error(
        "Could not open file " + path.string() + " for parsing");
  std::string contents(
      (std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
  return parse_matrix_value(contents);
}

inline std::vector<std::size_t> parse_indices(const std::string& s) {
  std::vector<std::string> tokens = utils::split_string(s);
  std::vector<std::size_t> indices(tokens.size());
  std::ranges::transform(tokens, indices.begin(), [](const std::string& t) {
    return static_cast<std::size_t>(std::stoull(t));
  });
  return indices;
}

inline std::vector<std::size_t> parse_indices_file(
    const std::filesystem::path& path) {
  std::ifstream file(path);
  if (!file.is_open())
    throw std::runtime_error(
        "Could not open file " + path.string() + " for parsing");
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

inline std::pair<std::string, std::vector<std::string>> parse_function(
    const std::string& s) {
  int paren_begin = s.find('(');
  int paren_end = s.find(')');

  // this is the case where f is used to denote f()
  if (paren_begin == std::string::npos || paren_end == std::string::npos)
    return {s, {}};

  // broken string
  if (paren_begin >= paren_end)
    throw std::runtime_error("Invalid functional form: " + s);

  std::string fn_name = s.substr(0, paren_begin);
  std::string params_str =
      s.substr(paren_begin + 1, paren_end - paren_begin - 1);
  std::vector<std::string> params = utils::split_string(params_str, ",");

  boost::to_lower(fn_name);

  return {fn_name, params};
}

inline amsim::Distribution parse_distribution(
    const std::string& s, bool is_probability = false) {
  auto [dist_name, dist_params_str] = parse_function(s);
  std::vector<double> dist_params(dist_params_str.size());

  if (!dist_params_str.empty())
    std::ranges::transform(
        dist_params_str, dist_params.begin(), [](const std::string& s) {
          return std::stod(s);
        });

  return amsim::str_to_distribution(dist_name, dist_params, is_probability);
}

inline amsim::WeightFunction parse_weight_function(const std::string& s) {
  auto [weight_name, weight_params_str] = parse_function(s);

  if (weight_name == "uniform") return amsim::Uniform();
  if (weight_name == "logistic") {
    Eigen::VectorXd logistic_weights(weight_params_str.size());
    std::ranges::transform(
        weight_params_str, logistic_weights.begin(), [](const std::string& s) {
          return std::stod(s);
        });
    return amsim::Logistic(logistic_weights);
  }

  throw std::runtime_error("Unrecognise weight function " + s);
}

}  // namespace amsim
