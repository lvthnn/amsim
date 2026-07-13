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

#include <amsim/core/log.h>

#include <Eigen/Dense>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace amsim {

struct GenomeParams {
  std::size_t n_ind;
  std::size_t n_loc;
  Eigen::VectorXd v_maf;
  Eigen::VectorXd v_rec;
  Eigen::VectorXd v_mut;
};

struct PhenomeParams {
  std::size_t n_pheno;              ///< Number of phenotypes
  std::vector<std::string> names;   ///< Vector of phenotype names
  std::vector<std::size_t> n_locs;  ///< Vector of number of loci per phenotype
  std::unordered_map<std::string, std::size_t> pheno_ids; ///< Name-ID map

  std::vector<Eigen::VectorXd> pheno_effects;  ///< Matrix of effect vectors
  std::vector<std::vector<std::size_t>> pheno_loc;  ///< Phenotype causal loci

  Eigen::VectorXd h2_gen;  ///< Vector of genetic component variances
  Eigen::VectorXd h2_env;  ///< Vector of environmental component variances
  Eigen::VectorXd h2_vert;  ///< Vector of nurture component variances

  Eigen::MatrixXd gen_cor;  ///< Genetic component correlation matrix
  Eigen::MatrixXd env_cor;  ///< Environmental component correlation matrix

  Eigen::VectorXd rnur_pat;  ///< Paternal ratio in nurture effect
  Eigen::VectorXd rnur_env;  ///< Environmental ratio in nurture effect
  Eigen::VectorXd vert_pat;  ///< Paternal vertical transmission ratio
};

struct MatingParams {
  Eigen::MatrixXd mate_cor;
  double tol_inf;
  std::size_t max_itr;
  double temp_init;
  double temp_decay;
};

struct SimulationParams {
  std::size_t n_gens;
  std::uint64_t rng_seed;
  std::optional<std::uint64_t> post_init_seed;
  std::filesystem::path out_dir;
  std::size_t rep_id;

  LogLevel log_level;
  bool log_to_file = false;
};

struct Params {
  GenomeParams geno;
  PhenomeParams pheno;
  MatingParams mate;
  SimulationParams sim;
};
}  // namespace amsim
