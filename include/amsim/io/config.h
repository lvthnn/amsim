#pragma once

#include <amsim/init/setup.h>
#include <toml++/toml.h>

#include <fstream>
#include <numeric>
#include <optional>

namespace amsim {

template <typename V>
inline std::string vec_to_string(const V& vec) {
  return std::accumulate(
      vec.begin() + 1,
      vec.end(),
      std::to_string(vec[0]),
      [](const std::string& a, double b) {
        return a + "," + std::to_string(b);
      });
}

inline std::string to_string(
    const std::variant<double, File<Eigen::MatrixXd>, Distribution>& val) {
  if (std::holds_alternative<double>(val))
    return std::to_string(std::get<double>(val));

  if (std::holds_alternative<Distribution>(val)) {
    const auto& dist = std::get<Distribution>(val);
    return dist.name + "(" + vec_to_string(dist.params) + ")";
  }

  return "file(" + std::get<File<Eigen::MatrixXd>>(val).path.string() + ")";
}

inline toml::array to_toml_array(const Eigen::MatrixXd& mat) {
  toml::array rows;
  for (int i = 0; i < mat.rows(); ++i) {
    toml::array row;
    for (int j = 0; j < mat.cols(); ++j) row.push_back(mat(i, j));
    rows.push_back(row);
  }
  return rows;
}

inline void insert_matrix_or_file(
    const std::variant<File<Eigen::MatrixXd>, Eigen::MatrixXd>& variant,
    const std::string& name,
    toml::table& table) {
  if (std::holds_alternative<File<Eigen::MatrixXd>>(variant))
    table.insert(name, to_string(std::get<File<Eigen::MatrixXd>>(variant)));
  else
    table.insert(name, to_toml_array(std::get<Eigen::MatrixXd>(variant)));
}

inline void insert_list_or_file(
    const std::variant<
        std::vector<std::size_t>,
        File<std::vector<std::size_t>>>& variant,
    const std::string& name,
    toml::table& table) {
  if (std::holds_alternative<File<std::vector<std::size_t>>>(variant)) {
    table.insert(
        name,
        "file(" +
            std::get<File<std::vector<std::size_t>>>(variant).path.string() +
            ")");
  } else {
    toml::array arr;
    for (auto idx : std::get<std::vector<std::size_t>>(variant))
      arr.push_back(static_cast<int64_t>(idx));
    table.insert(name, arr);
  }
}

inline void write_config(
    const Simulation& simulation,
    std::size_t n_replicates,
    std::size_t n_threads,
    const std::filesystem::path& config_file) {
  auto config = toml::table{
      {"simulation", toml::table{}},
      {"genome", toml::table{}},
      {"phenotypes", toml::table{}},
      {"mating", toml::table{}},
      {"estimators", toml::table{}}};

  // GLOBAL SIMULATION OPTIONS

  auto& sim = *config["simulation"].as_table();
  sim.insert("n_individuals", static_cast<int64_t>(simulation.n_individuals));
  sim.insert("n_generations", static_cast<int64_t>(simulation.n_generations));
  sim.insert("n_replications", static_cast<int64_t>(n_replicates));
  sim.insert("n_threads", static_cast<int64_t>(n_threads));

  if (simulation.random_seed.has_value())
    sim.insert("random_seed", static_cast<int64_t>(simulation.random_seed.value()));

  sim.insert("output_dir", simulation.output_dir.string());

  if (simulation.output_name.has_value())
    sim.insert("output_name", simulation.output_name.value());

  sim.insert("log_level", LogLevel_to_string(simulation.log_level));
  sim.insert("log_to_output", !simulation.log_to_file);

  // GENOME OPTIONS
  auto& genome = *config["genome"].as_table();
  genome.insert("n_loci", static_cast<int64_t>(simulation.genome.n_loci));
  genome.insert("locus_maf", to_string(simulation.genome.v_maf));
  genome.insert("locus_rec", to_string(simulation.genome.v_rec));
  genome.insert("locus_mut", to_string(simulation.genome.v_mut));

  // PHENOME OPTIONS
  auto& phenome = *config["phenotypes"].as_table();

  if (simulation.genetic_component_cor.has_value())
    insert_matrix_or_file(
        simulation.genetic_component_cor.value(), "genetic_cor", phenome);

  if (simulation.environmental_component_cor.has_value())
    insert_matrix_or_file(
        simulation.environmental_component_cor.value(),
        "environmental_cor",
        phenome);
  auto phenotypes = toml::array{};
  for (const auto& phenotype : simulation.phenotypes) {
    toml::table pheno{};
    pheno.insert("name", phenotype.name);
    pheno.insert("n_loci", static_cast<int64_t>(phenotype.n_causal_loci));
    if (phenotype.effects.has_value())
      pheno.insert("effects", to_string(phenotype.effects.value()));
    if (phenotype.causal_loci.has_value())
      insert_list_or_file(phenotype.causal_loci.value(), "causal_loci", pheno);
    pheno.insert("var_genetic", phenotype.var_genetic);
    pheno.insert("var_environmental", phenotype.var_environmental);
    pheno.insert("var_vertical", phenotype.var_vertical);
    phenotypes.push_back(pheno);
  }

  // mating — if it is assortative, need to log some extra information
  if (simulation.mating.type == "assortative") {
    auto& mating = *config["mating"].as_table();
    if (simulation.mating.mate_cor.has_value())
      insert_matrix_or_file(
          simulation.mating.mate_cor.value(), "mate_cor", mating);
    mating.insert("tol_inf", simulation.mating.tolerance);
    mating.insert("max_itr", static_cast<int64_t>(simulation.mating.max_iterations));
    mating.insert("temp_init", simulation.mating.initial_temperature);
    mating.insert("temp_decay", simulation.mating.temperature_decay);
  }

  // estimators
  auto& estimators = *config["estimators"].as_table();

  toml::array estimator_names;
  for (const auto& estimator : simulation.estimators)
    estimator_names.push_back(estimator.name);

  std::ofstream out(config_file);

  if (!out.is_open())
    throw std::runtime_error("Could not open configuration file for writing");

  out << config;
}

inline void read_config(const std::filesystem::path& config_file) {}

}  // namespace amsim
