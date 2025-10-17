#ifndef AMSIM_SIMULATION_PROGRESS_H
#define AMSIM_SIMULATION_PROGRESS_H

namespace amsim {

  enum class SimulationStatus {
    NEW        = 0,
    SIMULATION = 1,
    GENOME     = 2,
    PHENOME    = 3,
    MATING     = 4,
    READY      = 5,
    FINISHED   = 6
  };

  inline SimulationStatus operator++(SimulationStatus &status, int) {
    return SimulationStatus(int(status) + 1);
  }
}

#endif
