#ifndef AMSIMCPP_SIMULATION_PROGRESS_H
#define AMSIMCPP_SIMULATION_PROGRESS_H

namespace amsim {

enum class SimulationStatus {
  NEW = 0,
  GENOME = 1,
  PHENOME = 2,
  MATING = 3,
  READY = 4,
  FINISHED = 5
};

inline SimulationStatus operator++(SimulationStatus& status, int) {
  status = SimulationStatus(int(status) + 1);
  return status;
}

}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_PROGRESS_H
