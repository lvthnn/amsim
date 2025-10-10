#ifndef AMSIMCPP_SIMULATION_H
#define AMSIMCPP_SIMULATION_H

#include <amsim/genome.h>
#include <amsim/phenobuf.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/metric.h>

namespace amsim {

  class Simulation {
  public:
    Simulation(Genome genome_, PhenoArch);

  private:
    Genome genome_;
    PhenoArch arch_;
    PhenoBuf buf_;
    MatingModel model_;
    std::vector<Phenotype> phenotypes_;
  };

}

#endif //AMSIMCPP_SIMULATION_H