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

#include <amsim/sample/proband.h>
#include <amsim/estimate/sample.h>
#include <amsim/estimate/gwas.h>
#include <amsim/estimate/haseman_elston.h>
#include <amsim/estimate/external.h>
#include <amsim/estimate/greml.h>

namespace amsim {

template <Proband P>
inline SampleEstimator<P> make_sample_estimator(
    const SampleEstimatorDecl& decl) {
  if (decl.type == "external") {
    return SampleExternalEstimator<P>(
        decl.name,
        decl.exec.value(),
        decl.n_rows.value(),
        decl.n_cols.value(),
        decl.row_names,
        decl.col_names);
  }
  if (decl.type == "sample-mean") return SampleMeanEstimator<P>();
  if (decl.type == "sample-var") return SampleVarEstimator<P>();
  if (decl.type == "sample-cov") return SampleCovEstimator<P>();
  if (decl.type == "sample-mate-cor") return SampleMateCorEstimator<P>();
  if (decl.type == "haseman-elston")
    return SampleHasemanElstonEstimator<P>(decl.name);
  if (decl.type == "greml") return SampleGREMLEstimator<P>(decl.name);
  if (decl.type == "gwas") {
    std::size_t n_pcs = 0;
    double pval_threshold = 5e-8;

    if (decl.params.size() == 1)
      n_pcs = std::stoull(decl.params[0]);

    if (decl.params.size() == 2) {
      n_pcs = std::stoull(decl.params[0]);
      pval_threshold = std::stod(decl.params[1]);
    }

    return SampleGWASEstimator<P>(decl.name, n_pcs, pval_threshold);
  }
  throw std::runtime_error("Unknown sample estimator type: " + decl.type);
}

}  // namespace amsim
