#pragma once

#include <amsim/estimate/greml.h>
#include <amsim/estimate/gwas.h>
#include <amsim/estimate/haseman_elston.h>
#include <amsim/estimate/population.h>
#include <amsim/estimate/sample.h>
#include <amsim/estimate/external.h>

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
    if (auto it = decl.params.find("n_pcs"); it != decl.params.end())
      n_pcs = std::stoull(it->second);
    if (auto it = decl.params.find("pval_threshold"); it != decl.params.end())
      pval_threshold = std::stod(it->second);
    return SampleGWASEstimator<P>(decl.name, n_pcs, pval_threshold);
  }
  throw std::runtime_error("Unknown sample estimator type: " + decl.type);
}

}  // namespace amsim
