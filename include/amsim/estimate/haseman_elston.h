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
#include <amsim/sample/proband.h>

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <utility>

namespace amsim {

template <Proband P>
class HasemanElstonEstimator : public SampleEstimatorStrategy<P> {
 public:
  HasemanElstonEstimator(
      const Params& params,
      std::string sample_name,
      std::filesystem::path sample_dir,
      std::string name = "haseman-elston")
      : SampleEstimatorStrategy<P>(
            std::move(name),
            std::move(sample_name),
            sample_dir,
            params.pheno.names,
            std::vector<std::string>{"V(G)/Vp"},
            params.pheno.n_pheno,
            1),
        n_pheno_(params.pheno.n_pheno) {
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
      auto pfix = this->sample_dir_ / std::format("he_{}", p);
      utils::system_throttled(std::format(
          "gcta64 --grm {} --pheno {} --mpheno {} --HEreg "
          "--out {} --thread-num 1",
          (this->sample_dir_ / "grm").string(),
          (this->sample_dir_ / "data.pheno").string(),
          p + 1,
          pfix.string()));

      this->data_(p, 0) = parseHEreg(pfix.string() + ".HEreg");
    }
  }

 private:
  std::size_t n_pheno_;

  double parseHEreg(const std::string& path) {
    std::ifstream f(path);
    if (!f) return std::numeric_limits<double>::quiet_NaN();

    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
      std::istringstream ss(line);
      std::string key;
      double val;
      ss >> key >> val;
      if (key == "V(G)/Vp") return val;
    }
    return std::numeric_limits<double>::quiet_NaN();
  }
};

template <Proband P>
inline SampleEstimator<P> SampleHasemanElstonEstimator(
    std::string name = "haseman-elston") {
  return SampleEstimator<P>{
      .name = name,
      .fn = [name = std::move(name)](
                const Params& params,
                std::size_t /*n_probands*/,
                const std::filesystem::path& sample_dir) {
        return std::make_unique<HasemanElstonEstimator<P>>(
            params, sample_dir.filename().string(), sample_dir, name);
      }};
}

}  // namespace amsim
