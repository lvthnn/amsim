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
#include <amsim/core/utils.h>
#include <amsim/sample/proband.h>
#include <amsim/sample/weight.h>

#include <Eigen/Dense>
#include <algorithm>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <filesystem>
#include <fstream>
#include <vector>

namespace amsim {

inline std::string parse_exception_str(
    const std::string& s, const std::optional<std::string>& flag) {
  return std::format(

      "Failed to parse '{}' {}",
      s,
      flag.has_value() ? "(passed to " + flag.value() + ")" : "");
}

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

template <typename T>
inline T parse(const std::string& s);

template <typename T>
inline std::vector<T> parse_each(const std::vector<std::string>& ss) {
  std::vector<T> result(ss.size());
  for (std::size_t el = 0; el < ss.size(); ++el) result[el] = parse<T>(ss[el]);
  return result;
}

template <typename T>
inline std::vector<T> parse_vector(const std::string& s) {
  return parse_each<T>(utils::split_string(s));
}

template <typename T>
inline T parse(const std::string& s, const std::optional<std::string>& flag) {
  try {
    return parse<T>(s);
  } catch (const std::exception& e) {
    throw std::runtime_error(parse_exception_str(s, flag));
  }
}

template <>
inline std::size_t parse(const std::string& s) {
  return std::stoull(s);
}

template <>
inline std::uint64_t parse(const std::string& s) {
  return std::stoull(s);
}

template <>
inline double parse(const std::string& s) {
  return std::stod(s);
}

template <>
inline std::vector<std::size_t> parse(const std::string& s) {
  return parse_vector<std::size_t>(s);
}

template <>
inline std::vector<double> parse(const std::string& s) {
  return parse_vector<double>(s);
}

template <>
inline Eigen::VectorXd parse(const std::string& s) {
  std::vector<double> v = parse<std::vector<double>>(s);
  return Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
}

template <>
inline Eigen::MatrixXd parse(const std::string& s) {
  std::string norm = s;
  std::ranges::replace(norm, ';', '\n');
  std::ranges::replace(norm, '\t', ' ');
  std::ranges::replace(norm, ',', ' ');

  std::vector<std::string> row_strs = utils::split_string(norm, "\n");
  std::vector<Eigen::VectorXd> rows(row_strs.size());
  for (std::size_t r = 0; r < row_strs.size(); ++r)
    rows[r] = parse<Eigen::VectorXd>(row_strs[r]);

  std::size_t n_rows = rows.size();
  std::size_t n_cols = rows[0].size();
  auto bad = std::ranges::find_if(
      rows, [n_cols](const auto& r) { return r.size() != n_cols; });

  if (bad != rows.end())
    throw std::runtime_error(
        std::format(
            "Matrix row {} has {} values, expected {}",
            std::distance(rows.begin(), bad),
            bad->size(),
            n_cols));

  Eigen::MatrixXd matrix(n_rows, n_cols);

  for (std::size_t r = 0; r < rows.size(); ++r) {
    if (static_cast<std::size_t>(rows[r].size()) != n_cols)
      throw std::runtime_error(
          std::format(
              "Matrix row {} has {} values, expected {}",
              r,
              rows[r].size(),
              n_cols));
    matrix.row(r) = rows[r];
  }

  return (rows.size() == 1) ? matrix.transpose() : matrix;
}

template <>
inline Distribution parse(const std::string& s) {
  auto [dist_name, dist_params_str] = parse_function(s);
  std::vector<double> dist_params(dist_params_str.size());

  if (!dist_params_str.empty())
    dist_params = parse_each<double>(dist_params_str);

  return amsim::str_to_distribution(dist_name, dist_params);
}

template <>
inline WeightFunction parse(const std::string& s) {
  auto [weight_name, params_str] = parse_function(s);

  if (weight_name == "uniform") return amsim::Uniform();
  if (weight_name == "logistic") {
    auto params = parse_each<double>(params_str);
    return amsim::Logistic(
        Eigen::Map<Eigen::VectorXd>(params.data(), params.size()));
  }

  throw std::runtime_error("Unrecognised weight function " + s);
}

template <typename T>
inline std::vector<T> parse_vector(
    const std::string& s, const std::optional<std::string>& flag) {
  try {
    return parse_vector<T>(s);
  } catch (const std::exception& e) {
    throw std::runtime_error(parse_exception_str(s, flag));
  }
}

template <typename T>
inline T parse_file(const std::filesystem::path& path) {
  std::ifstream file(path);
  if (!file.is_open())
    throw std::runtime_error(
        "Could not open file " + path.string() + " for parsing");
  std::string contents(
      (std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
  return parse<T>(contents);
}

template <typename T>
struct File {
  std::filesystem::path path;
  T load() const { return parse_file<T>(path); }
};

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
    for (std::size_t c = 0; c < n_cols; ++c) result(r, c) = matrix[r][c];
  return result;
}

}  // namespace amsim
