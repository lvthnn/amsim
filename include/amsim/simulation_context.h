#ifndef AMSIMCPP_SIMULATION_CONTEXT_H
#define AMSIMCPP_SIMULATION_CONTEXT_H

#include <amsim/genome.h>
#include <amsim/mating.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/phenotype.h>

#include <functional>

namespace amsim {

struct SimulationContext {
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

  // general data
  const std::size_t n_ind;
  const std::size_t n_sex;
  const std::size_t n_loc;
  const std::size_t n_pheno;
  const std::vector<std::string> pheno_names;

  // simulation components
  Genome& genome;
  PhenoArch& arch;
  PhenoBuf& buf;
  PhenotypeList& phenotypes;
  AssortativeModel& model;
};

using MetricFunc = std::function<std::vector<double>(const SimulationContext&)>;

}  // namespace amsim

#endif  // AMSIM_SIMULATION_CONTEXT_H
