#ifndef AMSIMCPP_SIMULATION_H
#define AMSIMCPP_SIMULATION_H

#include <string>
#include <fstream>

#include <amsim/genome.h>
#include <amsim/phenobuf.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/mating.h>
#include <amsim/metric.h>

#include <amsim/simulation_status.h>

namespace amsim {

  class Simulation {
  public:
    Simulation(std::size_t n_gen, std::string out_dir, Genome &genome,
               PhenoArch &arch, PhenoBuf &buf, AssortativeModel &model,
               std::vector<Metric> &metrics);
    void run();

  private:
    const std::size_t n_gen_;
    const std::string out_dir_;
    SimulationStatus status_;

    Genome            genome_;
    PhenoArch         arch_;
    PhenoBuf          buf_;
    PhenotypeList     phenotypes_;
    AssortativeModel  model_;

    std::vector<Metric> metrics_;
    std::vector<std::ofstream> streams_;

    void stream_();
  };

}

#endif //AMSIMCPP_SIMULATION_H
