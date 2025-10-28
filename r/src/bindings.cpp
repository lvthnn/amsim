// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>

#include <string>

#include <amsim/component_type.h>
#include <amsim/mating_type.h>
#include <amsim/metricspec.h>
#include <amsim/simulation_builder.h>

//-----------------------------------------------------------------------------
// enum bindings
//-----------------------------------------------------------------------------

amsim::ComponentType _ComponentType(std::string component_type) {
  if (component_type == "genetic")
    return amsim::ComponentType::GENETIC;
  if (component_type == "environmental")
    return amsim::ComponentType::ENVIRONMENTAL;
  if (component_type == "vertical")
    return amsim::ComponentType::VERTICAL;
  if (component_type == "total")
    return amsim::ComponentType::TOTAL;
  Rcpp::stop("Unknown component type '" + component_type + "'");
}

amsim::MatingType _MatingType(std::string mating_type) {
  if (mating_type == "random")
    return amsim::MatingType::RANDOM;
  if (mating_type == "assortative")
    return amsim::MatingType::ASSORTATIVE;
  Rcpp::stop("Unknown mating type '" + mating_type + "'");
}

//-----------------------------------------------------------------------------
// expose metric specs as external pointers
//-----------------------------------------------------------------------------

SEXP _def_metric(amsim::MetricSpec spec) {
  return Rcpp::XPtr<amsim::MetricSpec>(new amsim::MetricSpec(spec));
}

// [[Rcpp::export(".pheno_h2")]]
SEXP _pheno_h2() {
  return _def_metric(amsim::pheno_h2());
}

// [[Rcpp::export(".pheno_comp_mean")]]
SEXP _pheno_comp_mean(std::string component_type) {
  return _def_metric(amsim::pheno_comp_mean(_ComponentType(component_type)));
}

// [[Rcpp::export(".pheno_comp_var")]]
SEXP _pheno_comp_var(std::string component_type) {
  return _def_metric(amsim::pheno_comp_var(_ComponentType(component_type)));
}

// [[Rcpp::export(".pheno_comp_cor")]]
SEXP _pheno_comp_cor(std::string component_type) {
  return _def_metric(amsim::pheno_comp_cor(_ComponentType(component_type)));
}

// [[Rcpp::export(".pheno_comp_xcor")]]
SEXP _pheno_comp_xcor(std::string component_type) {
  return _def_metric(amsim::pheno_comp_xcor(_ComponentType(component_type)));
}

// [[Rcpp::export(".pheno_latent_h2")]]
SEXP _pheno_latent_h2() {
  return _def_metric(amsim::pheno_latent_h2());
}

// [[Rcpp::export(".pheno_latent_comp_mean")]]
SEXP _pheno_latent_comp_mean(std::string component_type) {
  return _def_metric(amsim::pheno_latent_comp_mean(_ComponentType(component_type)));
}

// [[Rcpp::export(".pheno_latent_comp_var")]]
SEXP _pheno_latent_comp_var(std::string component_type) {
  return _def_metric(amsim::pheno_latent_comp_var(_ComponentType(component_type)));
}

// [[Rcpp::export(".pheno_latent_comp_cor")]]
SEXP _pheno_latent_comp_cor(std::string component_type) {
  return _def_metric(amsim::pheno_latent_comp_cor(_ComponentType(component_type)));
}

// [[Rcpp::export(".pheno_latent_comp_xcor")]]
SEXP _pheno_latent_comp_xcor(std::string component_type) {
  return _def_metric(amsim::pheno_latent_comp_xcor(_ComponentType(component_type)));
}

//-----------------------------------------------------------------------------
// expose simulation components to R bindings
//-----------------------------------------------------------------------------

// [[Rcpp::export(".builder_new")]]
SEXP _builder_new() {
  amsim::SimulationBuilder* builder = new amsim::SimulationBuilder();
  return Rcpp::XPtr<amsim::SimulationBuilder>(builder);
}

void _builder_simulation() {

}
