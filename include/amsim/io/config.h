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

#include <amsim/init/setup.h>
#include <toml++/toml.h>

#include <fstream>
#include <numeric>

namespace amsim {

template <typename V>
inline std::string vec_to_string(const V& vec) {
  return std::accumulate(
      vec.begin() + 1,
      vec.end(),
      std::format("{:g}", vec[0]),
      [](const std::string& a, double b) {
        return a + "," + std::format("{:g}", b);
      });
}

inline std::string to_string(
    const std::variant<double, File<Eigen::MatrixXd>, Distribution>& val) {
  if (std::holds_alternative<double>(val))
    return std::format("{:g}", std::get<double>(val));

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

inline toml::array to_toml_array(const std::vector<std::string>& vec) {
  toml::array arr;
  for (const auto& s : vec) arr.push_back(s);
  return arr;
}

inline void insert_const_file_or_dist(
    toml::table& table,
    const std::string& name,
    const std::variant<double, File<Eigen::MatrixXd>, Distribution>& variant) {
  if (std::holds_alternative<double>(variant)) {
    table.insert(name, std::get<double>(variant));
    return;
  }
  if (std::holds_alternative<File<Eigen::MatrixXd>>(variant)) {
    table.insert(
        name,
        "file(" + std::get<File<Eigen::MatrixXd>>(variant).path.string() + ")");
    return;
  }

  Distribution dist = std::get<Distribution>(variant);
  table.insert(name, dist.name + "(" + vec_to_string(dist.params) + ")");
}

inline void insert_matrix_or_file(
    toml::table& table,
    const std::string& name,
    const std::variant<File<Eigen::MatrixXd>, Eigen::MatrixXd>& variant) {
  if (std::holds_alternative<File<Eigen::MatrixXd>>(variant))
    table.insert(
        name,
        "file(" + std::get<File<Eigen::MatrixXd>>(variant).path.string() + ")");
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
    const std::vector<SampleDecl>& sample_decl,
    const std::vector<SampleEstimatorDecl>& sample_estimator_decl,
    const std::filesystem::path& config_file) {
  auto config = toml::table{
      {"simulation", toml::table{}},
      {"genome", toml::table{}},
      {"phenotypes", toml::table{}},
      {"mating", toml::table{}},
      {"estimators", toml::table{}},
      {"samples", toml::table{}}};

  // GLOBAL SIMULATION OPTIONS

  auto& sim = *config["simulation"].as_table();
  sim.insert("n_individuals", static_cast<int64_t>(simulation.n_individuals));
  sim.insert("n_generations", static_cast<int64_t>(simulation.n_generations));
  sim.insert("n_replications", static_cast<int64_t>(n_replicates));
  sim.insert("n_threads", static_cast<int64_t>(n_threads));

  if (simulation.random_seed.has_value())
    sim.insert(
        "random_seed", static_cast<int64_t>(simulation.random_seed.value()));

  sim.insert("output_dir", simulation.output_dir.string());

  if (simulation.output_name.has_value())
    sim.insert("output_name", simulation.output_name.value());

  sim.insert("log_level", LogLevel_to_string(simulation.log_level));
  sim.insert("log_to_output", !simulation.log_to_file);

  // GENOME OPTIONS
  auto& genome = *config["genome"].as_table();
  genome.insert("n_loci", static_cast<int64_t>(simulation.genome.n_loci));
  insert_const_file_or_dist(genome, "locus_maf", simulation.genome.v_maf);
  insert_const_file_or_dist(genome, "locus_rec", simulation.genome.v_rec);
  insert_const_file_or_dist(genome, "locus_mut", simulation.genome.v_mut);

  // PHENOME OPTIONS
  auto& phenotype_sec = *config["phenotypes"].as_table();

  if (simulation.genetic_component_cor.has_value())
    insert_matrix_or_file(
        phenotype_sec, "genetic_cor", simulation.genetic_component_cor.value());

  if (simulation.environmental_component_cor.has_value())
    insert_matrix_or_file(
        phenotype_sec,
        "environmental_cor",
        simulation.environmental_component_cor.value());

  auto phenotypes = toml::array{};
  for (const auto& phenotype : simulation.phenotypes) {
    toml::table phenotype_tbl{};
    phenotype_tbl.insert("name", phenotype.name);
    phenotype_tbl.insert(
        "n_loci", static_cast<int64_t>(phenotype.n_causal_loci));
    if (phenotype.effects.has_value())
      phenotype_tbl.insert("effects", to_string(phenotype.effects.value()));
    if (phenotype.causal_loci.has_value())
      insert_list_or_file(
          phenotype.causal_loci.value(), "causal_loci", phenotype_tbl);
    phenotype_tbl.insert("var_genetic", phenotype.var_genetic);
    phenotype_tbl.insert("var_environmental", phenotype.var_environmental);
    phenotype_tbl.insert("var_vertical", phenotype.var_vertical);
    phenotypes.push_back(phenotype_tbl);
  }
  phenotype_sec.insert("phenotype", phenotypes);

  // mating — if it is assortative, need to log some extra information
  if (simulation.mating.type == "assortative") {
    auto& mating = *config["mating"].as_table();
    if (simulation.mating.mate_cor.has_value())
      insert_matrix_or_file(
          mating, "mate_cor", simulation.mating.mate_cor.value());
    mating.insert("tol_inf", simulation.mating.tolerance);
    mating.insert(
        "max_itr", static_cast<int64_t>(simulation.mating.max_iterations));
    mating.insert("temp_init", simulation.mating.initial_temperature);
    mating.insert("temp_decay", simulation.mating.temperature_decay);
  }

  // estimators
  auto& estimators = *config["estimators"].as_table();

  toml::array estimator_names;
  for (const auto& estimator : simulation.estimators)
    estimator_names.push_back(estimator.name);

  estimators.insert("population", estimator_names);

  // samples
  if (!sample_decl.empty()) {
    auto& sample_section = *config["samples"].as_table();

    toml::array samples{};
    for (const auto& sample : sample_decl) {
      toml::table sample_tbl{};
      sample_tbl.insert("name", sample.name);
      sample_tbl.insert("proband_type", sample.proband_type);
      sample_tbl.insert("n_probands", static_cast<int64_t>(sample.n_probands));
      if (sample.on.has_value())
        sample_tbl.insert("on", to_toml_array(sample.on.value()));
      if (sample.of.has_value())
        sample_tbl.insert("of", to_toml_array(sample.of.value()));
      if (sample.weight_function.has_value())
        sample_tbl.insert("weight_function", sample.weight_function.value());
      if (sample.agg.has_value())
        sample_tbl.insert("agg", sample.agg.value());
      if (!sample.estimators.empty())
        sample_tbl.insert("estimators", to_toml_array(sample.estimators));
      samples.push_back(sample_tbl);
    }
    sample_section.insert("sample", samples);

    toml::array sample_estimators{};
    for (const auto& sample_estimator : sample_estimator_decl) {
      toml::table sample_estimator_tbl{};
      sample_estimator_tbl.insert("name", sample_estimator.name);
      sample_estimator_tbl.insert("type", sample_estimator.type);
      if (sample_estimator.exec.has_value()) {
        sample_estimator_tbl.insert("exec", sample_estimator.exec.value());
        if (sample_estimator.n_rows.has_value())
          sample_estimator_tbl.insert(
              "n_rows", static_cast<int64_t>(sample_estimator.n_rows.value()));
        if (sample_estimator.n_cols.has_value())
          sample_estimator_tbl.insert(
              "n_cols", static_cast<int64_t>(sample_estimator.n_cols.value()));
        if (sample_estimator.row_names.has_value())
          sample_estimator_tbl.insert(
              "row_names", to_toml_array(sample_estimator.row_names.value()));
        if (sample_estimator.col_names.has_value())
          sample_estimator_tbl.insert(
              "col_names", to_toml_array(sample_estimator.col_names.value()));
      }
      sample_estimators.push_back(sample_estimator_tbl);
    }
    sample_section.insert("sample_estimators", sample_estimators);
  }

  std::ofstream out(config_file);

  if (!out.is_open())
    throw std::runtime_error("Could not open configuration file for writing");

  out << toml::
             toml_formatter{config, toml::format_flags::relaxed_float_precision}
      << std::endl;
}

inline void read_config(const std::filesystem::path& config_file) {}

}  // namespace amsim
