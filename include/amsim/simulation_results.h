#ifndef AMSIMCPP_SIMULATION_RESULTS_H
#define AMSIMCPP_SIMULATION_RESULTS_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>

namespace amsim {

struct CellStats {
  CellStats(std::vector<double> stream);

  const double n_elem;
  const double mean;
  const double median;
  const double std;
  const double sem;
  const double lci;
  const double uci;
  const double q025;
  const double q975;
};

using KeyMap = std::unordered_map<std::string, int>;
using ResultsTable = std::vector<std::vector<std::vector<CellStats>>>;

class SimulationResults {
 public:
  SimulationResults(
      std::filesystem::path out_dir,
      std::optional<std::size_t> n_replicates = std::nullopt,
      std::optional<std::set<std::string>> metric_names = std::nullopt);

  void load();
  void summarise();
  void save(std::filesystem::path out_path);
  void print(std::string metric_name);

 private:
  std::filesystem::path out_dir_;
  std::size_t n_replicates_ = 0;
  std::size_t n_metrics_;

  std::set<std::string> metric_names_;
  std::vector<std::vector<std::string>> label_names_;
  std::vector<std::filesystem::path> rep_dirs_;

  KeyMap metric_map_;
  std::vector<KeyMap> label_map_;
  ResultsTable table_;

  void summarise_metric_(std::string metric_name);
};

}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_RESULTS_H
