#ifndef AMSIMCPP_SIMULATION_CONTEXT_H
#define AMSIMCPP_SIMULATION_CONTEXT_H

#pragma once
#include <functional>

#include <amsim/genome.h>
#include <amsim/phenotype.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/mating.h>

namespace amsim {

  struct SimulationContext {
    SimulationContext(Genome &genome_, PhenoArch &arch_, PhenoBuf &buf_,
                      PhenotypeList &phenotypes_, AssortativeModel &model_)
      : genome(genome_),
        arch(arch_),
        buf(buf_),
        phenotypes(phenotypes_),
        model(model_) {};

    Genome           &genome;
    PhenoArch        &arch;
    PhenoBuf         &buf;
    PhenotypeList    &phenotypes;
    AssortativeModel &model;
  };

  using MetricFunc = std::function<std::vector<double>(const SimulationContext&)>;
}

#endif // AMSIM_SIMULATION_CONTEXT_H
