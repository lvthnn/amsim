#ifndef AMSIMCPP_SIMULATION_H
#define AMSIMCPP_SIMULATION_H

#include <string>

#include <amsim/genome.h>
#include <amsim/phenobuf.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/mating.h>
#include <amsim/metric.h>

namespace amsim {

  class Simulation {
  public:
    Simulation(Genome &genome, PhenoArch &arch, PhenoBuf &buf,
               AssortativeModel &model);

    void run();

  private:
    Genome            genome_;
    PhenoArch         arch_;
    PhenoBuf          buf_;
    PhenotypeList     pheno_;
    AssortativeModel  model_;

    const std::size_t n_itr_;
    const std::string out_dir_;
  };

}

#endif //AMSIMCPP_SIMULATION_H
