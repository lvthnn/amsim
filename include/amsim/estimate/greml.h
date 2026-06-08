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

#include <amsim/core/params.h>
#include <amsim/core/utils.h>
#include <amsim/estimate/sample.h>
#include <amsim/io/writer.h>
#include <amsim/sample/proband.h>

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <utility>

namespace amsim {

template <Proband P>
class GREMLEstimator : public SampleEstimatorStrategy<P> {
 public:
  GREMLEstimator(
      const Params& params,
      std::string sample_name,
      std::filesystem::path sample_dir,
      std::string name = "greml")
      : SampleEstimatorStrategy<P>(
            std::move(name),
            std::move(sample_name),
            sample_dir,
            params.pheno.names,
            {"V(G)", "V(E)", "V(G)/[V(G) + V(E)]"},
            params.pheno.n_pheno,
            3),
        n_pheno_(params.pheno.n_pheno),
        pheno_names_(params.pheno.names) {
    utils::check_gcta64();
  }

  void compute() override {
    if (!std::filesystem::exists(this->sample_dir_ / "grm.grm.bin")) {
      utils::system_throttled(std::format(
          "gcta64 --bfile {} --make-grm --out {} "
          "--thread-num 1",
          (this->sample_dir_ / "data").string(),
          (this->sample_dir_ / "grm").string()));
    }

    for (std::size_t p = 0; p < n_pheno_; ++p) {
      auto pfix = this->sample_dir_ / std::format("greml_{}", p);
      utils::system_throttled(std::format(
          "gcta64 --grm {} --pheno {} --mpheno {} --reml "
          "--out {} --thread-num 1",
          (this->sample_dir_ / "grm").string(),
          (this->sample_dir_ / "data.pheno").string(),
          p + 1,
          pfix.string()));

      this->data_.row(p) = parseHsq(pfix.string() + ".hsq");
    }
  }

 private:
  std::size_t n_pheno_;
  std::vector<std::string> pheno_names_;

  Eigen::RowVector3d parseHsq(const std::string& path) {
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    Eigen::RowVector3d result{NaN, NaN, NaN};

    std::ifstream f(path);
    if (!f) return result;

    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
      std::istringstream ss(line);
      std::string key;
      double val;
      ss >> key >> val;
      if (key == "V(G)")
        result(0) = val;
      else if (key == "V(e)")
        result(1) = val;
      else if (key == "V(G)/Vp")
        result(2) = val;
    }
    return result;
  }
};

template <Proband P>
inline SampleEstimator<P> SampleGREMLEstimator(std::string name = "greml") {
  return SampleEstimator<P>{
      .name = name,
      .fn = [name = std::move(name)](
                const Params& params,
                std::size_t /*n_probands*/,
                const std::filesystem::path& sample_dir) {
        return std::make_unique<GREMLEstimator<P>>(
            params, sample_dir.filename().string(), sample_dir, name);
      }};
}

}  // namespace amsim
