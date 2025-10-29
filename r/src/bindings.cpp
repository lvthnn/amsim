// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <optional>

#include <amsim/component_type.h>
#include <amsim/mating_type.h>
#include <amsim/metricspec.h>
#include <amsim/simulation_builder.h>

//-----------------------------------------------------------------------------
// convenience functions 
//-----------------------------------------------------------------------------

template<typename T>
std::optional<T> _optional(SEXP x) {
  if (x == R_NilValue) return std::nullopt;
  return Rcpp::as<T>(x);
}

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

std::vector<amsim::MetricSpec> _MetricSpecs(Rcpp::List metrics) {
  std::size_t n_specs = metrics.length();
  std::vector<amsim::MetricSpec> metric_specs;
  metric_specs.reserve(n_specs);
  for (std::size_t spec = 0; spec < n_specs; spec++) {
    Rcpp::XPtr<amsim::MetricSpec> spec_cur(metrics[spec]);
    metric_specs.push_back(*spec_cur);
  }
  return metric_specs;
}

//-----------------------------------------------------------------------------
// expose metric specs as external pointers
//-----------------------------------------------------------------------------

SEXP _def_metric(amsim::MetricSpec spec) {
  return Rcpp::XPtr<amsim::MetricSpec>(new amsim::MetricSpec(spec), true);
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
  return Rcpp::XPtr<amsim::SimulationBuilder>(builder, true);
}

// [[Rcpp::export(".builder_simulation")]]
void _builder_simulation(
    SEXP sexp_builder,
    std::size_t n_generations,
    std::size_t n_individuals,
    std::string output_dir,
    std::uint64_t random_seed) {
  Rcpp::XPtr<amsim::SimulationBuilder> builder(sexp_builder); 
  builder->simulation(n_generations, n_individuals, output_dir, random_seed);
}

// [[Rcpp::export(".builder_genome")]]
void _builder_genome(
    SEXP sexp_builder,
    std::size_t n_loci,
    std::vector<double> locus_mafs,
    std::vector<double> locus_recombination,
    std::vector<double> locus_mutation) {
  Rcpp::XPtr<amsim::SimulationBuilder> builder(sexp_builder);
  builder->genome(n_loci, locus_mafs, locus_recombination, locus_mutation);
}

// [[Rcpp::export(".builder_phenome")]]
void _builder_phenome(
    SEXP sexp_builder,
    std::size_t n_phenotypes,
    std::vector<std::string> names,
    std::vector<std::size_t> loci,
    std::vector<double> h2_genetic,
    std::vector<double> h2_environmental,
    std::vector<double> h2_vertical,
    std::vector<double> genetic_cor,
    std::vector<double> environmental_cor) {
  Rcpp::XPtr<amsim::SimulationBuilder> builder(sexp_builder);
  builder->phenome(
      n_phenotypes,
      names,
      loci,
      h2_genetic,
      h2_environmental,
      h2_vertical,
      genetic_cor,
      environmental_cor);
}

// [[Rcpp::export(".builder_mating")]]
void _builder_mating(
    SEXP sexp_builder,
    std::string mating_type,
    SEXP n_iterations,
    SEXP temp_init,
    SEXP temp_decay,
    SEXP mate_cor) {
  Rcpp::XPtr<amsim::SimulationBuilder> builder(sexp_builder);
  builder->mating(
      _MatingType(mating_type),
      _optional<std::size_t>(n_iterations),
      _optional<double>(temp_init),
      _optional<double>(temp_decay),
      _optional<std::vector<double>>(mate_cor));
}

// [[Rcpp::export(".builder_metrics")]]
void _builder_metrics(SEXP sexp_builder, Rcpp::List metrics) {
  Rcpp::XPtr<amsim::SimulationBuilder> builder(sexp_builder);
  std::vector<amsim::MetricSpec> metric_specs(_MetricSpecs(metrics));
  builder->metrics(metric_specs);
}

// [[Rcpp::export(".builder_build")]]
SEXP _builder_build(SEXP sexp_builder) {
  Rcpp::XPtr<amsim::SimulationBuilder> builder(sexp_builder);
  amsim::Simulation* simulation = new amsim::Simulation(builder->build());
  return Rcpp::XPtr<amsim::Simulation>(simulation, true);
}

// [[Rcpp::export(".simulation_run")]]
void _simulation_run(SEXP sexp_simulation) {
  Rcpp::XPtr<amsim::Simulation> simulation(sexp_simulation);
  simulation->run();
}
