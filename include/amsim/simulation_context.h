#ifndef AMSIMCPP_SIMULATION_CONTEXT_H
#define AMSIMCPP_SIMULATION_CONTEXT_H

#include <amsim/genome.h>
#include <amsim/mating.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/phenotype.h>

#include <functional>

namespace amsim {

/// @brief Context holding references to simulation components
///
/// SimulationContext aggregates references to the genome, phenotypes, and
/// mating model, providing a single object to pass to metric functions.
struct SimulationContext {
  /// @brief Construct a SimulationContext
  ///
  /// @param genome_ Genome reference
  /// @param arch_ Phenotype architecture reference
  /// @param buf_ Phenotype buffer reference
  /// @param phenotypes_ Phenotype list reference
  /// @param model_ Mating model reference
  SimulationContext(
      Genome& genome_,
      PhenoArch& arch_,
      PhenoBuf& buf_,
      PhenotypeList& phenotypes_,
      AssortativeModel& model_)
      : n_ind(genome_.n_ind()),
        n_sex(genome_.n_ind() / 2),
        n_loc(genome_.n_loc()),
        n_pheno(phenotypes_.size()),
        pheno_names([&]() {
          std::vector<std::string> names;
          names.resize(phenotypes_.size());
          for (std::size_t pheno = 0; pheno < phenotypes_.size(); ++pheno)
            names[pheno] = phenotypes_[pheno].name();
          return names;
        }()),
        genome(genome_),
        arch(arch_),
        buf(buf_),
        phenotypes(phenotypes_),
        model(model_) {};

  const std::size_t n_ind;                     ///< Population size
  const std::size_t n_sex;                     ///< Individuals per sex
  const std::size_t n_loc;                     ///< Number of loci
  const std::size_t n_pheno;                   ///< Number of phenotypes
  const std::vector<std::string> pheno_names;  ///< Phenotype names

  Genome& genome;             ///< Genome reference
  PhenoArch& arch;            ///< Phenotype architecture reference
  PhenoBuf& buf;              ///< Phenotype buffer reference
  PhenotypeList& phenotypes;  ///< Phenotype list reference
  AssortativeModel& model;    ///< Mating model reference
};

/// @brief Function type for computing metrics from simulation context
using MetricFunc = std::function<std::vector<double>(const SimulationContext&)>;

}  // namespace amsim

#endif  // AMSIM_SIMULATION_CONTEXT_H
