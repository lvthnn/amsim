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

#include <amsim/estimate/estimator.h>
#include <amsim/estimate/genome_estimators.h>
#include <amsim/estimate/mating_estimators.h>
#include <amsim/estimate/pedigree_estimators.h>
#include <amsim/estimate/phenome_estimators.h>

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace amsim {

using PopulationEstimatorFactory =
    std::function<PopulationEstimator(const std::vector<std::string>&)>;

inline PopulationEstimator build_population_estimator(
    const std::string& name, const std::vector<std::string>& params) {
  static const std::unordered_map<std::string, PopulationEstimatorFactory>
      Registry = {
          {"genotype-mean",
           [](const std::vector<std::string>& /*s*/) {
             return PopulationGenotypeMean();
           }},
          {"genotype-var",
           [](const std::vector<std::string>& /*s*/) {
             return PopulationGenotypeVar();
           }},
          {"genotype-maf",
           [](const std::vector<std::string>& /*s*/) {
             return PopulationGenotypeMAF();
           }},
          {"genotype-cov",
           [](const std::vector<std::string>& /*s*/) {
             return PopulationGenotypeCov();
           }},
          {"genotype-cor",
           [](const std::vector<std::string>& /*s*/) {
             return PopulationGenotypeCor();
           }},
          {"heritability",
           [](const std::vector<std::string>& /*s*/) {
             return PopulationHeritability();
           }},
          {"pheno-mean",
           [](const std::vector<std::string>& s) {
             Component component = Component::Total;
             if (!s.empty()) component = Component_from_string(s[0]);
             return PopulationComponentMean(component);
           }},
          {"pheno-var",
           [](const std::vector<std::string>& s) {
             Component component = Component::Total;
             if (!s.empty()) component = Component_from_string(s[0]);
             return PopulationComponentMean(component);
           }},
          {"pheno-cor",
           [](const std::vector<std::string>& s) {
             Component component = Component::Total;
             if (!s.empty()) component = Component_from_string(s[0]);
             return PopulationComponentCor(component);
           }},
          {"pheno-cov",
           [](const std::vector<std::string>& s) {
             Component component = Component::Total;
             if (!s.empty()) component = Component_from_string(s[0]);
             return PopulationComponentCov(component);
           }},
          {"mate-cor",
           [](const std::vector<std::string>& s) {
             Component component = Component::Total;
             if (!s.empty()) component = Component_from_string(s[0]);
             return PopulationMateCor(component);
           }},
          {"cousin-cov",
           [](const std::vector<std::string>& s) {
             std::size_t degree = 1;
             Component component = Component::Total;
             if (!s.empty()) component = Component_from_string(s[0]);
             if (s.size() > 1) degree = std::stoull(s[1]);

             return PopulationCousinCov(degree, component);
           }},
          {"ancestor-cov", [](const std::vector<std::string>& s) {
             std::size_t degree = 1;
             Component component = Component::Total;
             if (!s.empty()) component = Component_from_string(s[0]);
             if (s.size() > 1) degree = std::stoull(s[1]);

             return PopulationAncestorCov(degree, component);
           }}};

  auto it = Registry.find(name);
  if (it == Registry.end())
    throw std::runtime_error("Could not find population estimator " + name);

  return it->second(params);
}

}  // namespace amsim
