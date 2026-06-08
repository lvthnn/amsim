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

#include <amsim/core/utils.h>
#include <amsim/estimate/sample.h>
#include <amsim/io/parse.h>
#include <amsim/sample/proband.h>

#include <utility>

namespace amsim {

template <Proband P>
class ExternalEstimator : public SampleEstimatorStrategy<P> {
 public:
  ExternalEstimator(
      std::string_view sample_name,
      const std::filesystem::path& sample_dir,
      std::string_view name,
      std::string_view exec,
      std::size_t n_rows,
      std::size_t n_cols = 1,
      std::optional<std::vector<std::string>> row_labels = std::nullopt,
      std::optional<std::vector<std::string>> col_labels = std::nullopt)
      : SampleEstimatorStrategy<P>(
            name,
            sample_name,
            sample_dir,
            std::move(row_labels.value_or(std::vector<std::string>{})),
            std::move(col_labels.value_or(std::vector<std::string>{})),
            n_rows,
            n_cols),
        exec_(exec),
        bfile_(sample_dir / "data"),
        out_path_(sample_dir / std::format("results_{}", name)) {}

  void compute() override {
    run_command();
    this->data_ = parse_matrix_file(out_path_);
  }

 private:
  std::string exec_;
  std::filesystem::path bfile_;
  std::filesystem::path out_path_;

  void run_command() {
    utils::system_throttled(
        std::format(
            "{} --bfile {} --out {}",
            exec_,
            bfile_.string(),
            out_path_.string()));
  }
};

template <Proband P>
inline SampleEstimator<P> SampleExternalEstimator(
    std::string name,
    std::string exec,
    std::size_t n_rows,
    std::size_t n_cols = 1,
    std::optional<std::vector<std::string>> row_labels = std::nullopt,
    std::optional<std::vector<std::string>> col_labels = std::nullopt) {
  return SampleEstimator<P>{
      .name = name,
      .fn = [name = std::move(name),
             exec = std::move(exec),
             n_rows,
             n_cols,
             row_labels = std::move(row_labels),
             col_labels = std::move(col_labels)](
                const Params& /*params*/,
                std::size_t /*n_probands*/,
                const std::filesystem::path& sample_dir) {
        return std::make_unique<ExternalEstimator<P>>(
            sample_dir.filename().string(),
            sample_dir,
            name,
            exec,
            n_rows,
            n_cols,
            row_labels,
            col_labels);
      }};
}

}  // namespace amsim
