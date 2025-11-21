#ifndef AMSIMCPP_SIMULATION_RESULTS_H
#define AMSIMCPP_SIMULATION_RESULTS_H

#include <amsim/stats.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace amsim {

constexpr std::size_t NUM_STATS = 8;
using KeyMap = std::unordered_map<std::string, std::size_t>;
using InvKeyMap = std::unordered_map<std::size_t, std::string>;

struct ResultsTable {
  ResultsTable() = default;

  void add(std::vector<double> stream) {
    std::size_t n_elem = stream.size();

    const double* ptr = stream.data();
    double mean = stats::mean(n_elem, ptr, 1);
    double sem = stats::sem(n_elem, ptr, 1);

    data[0].push_back(stats::mean(n_elem, ptr, 1));
    data[1].push_back(stats::quantile(0.5, n_elem, ptr, 1));
    data[2].push_back(stats::std(n_elem, ptr, 1, false));
    data[3].push_back(sem);
    data[4].push_back(mean + (1.96 * sem));
    data[5].push_back(mean - (1.96 * sem));
    data[6].push_back(stats::quantile(0.025, n_elem, ptr, 1));
    data[7].push_back(stats::quantile(0.975, n_elem, ptr, 1));
  }

  KeyMap label_map;
  InvKeyMap inv_label_map;
  std::array<std::vector<double>, NUM_STATS> data;
};

using ResultsIndex = std::vector<ResultsTable>;

class SimulationResults {
 public:
  explicit SimulationResults(
      std::filesystem::path out_dir,
      std::optional<std::size_t> n_replicates = std::nullopt,
      std::optional<std::vector<std::string>> metric_names = std::nullopt);

  ResultsTable operator()(const std::string& metric);

  void summarise();
  void save(
      std::optional<std::vector<std::string>> metrics = std::nullopt,
      std::optional<std::filesystem::path> out_dir = std::nullopt,
      bool overwrite = false);
  std::vector<std::string> metric_names() { return metric_names_; };

 private:
  // methods to use for constructor
  void inferReplicates(std::optional<std::size_t> n_replicates = std::nullopt);
  void inferMetrics(
      std::optional<std::vector<std::string>> metric_names = std::nullopt);
  static void encodeValues(
      const std::vector<std::string>& values,
      KeyMap& map_out,
      InvKeyMap& invmap_out);
  static void getLabels(
      std::vector<std::string>& labels, std::vector<std::ifstream>& streams);
  void summariseMetric(const std::string& metric);

  // inferred or supplied parameters
  std::filesystem::path out_dir_;

  // replicate fields
  std::size_t n_replicates_ = 0;
  std::vector<std::filesystem::path> rep_dirs_;

  // metric fields
  std::size_t n_metrics_;
  std::vector<std::string> metric_names_;
  std::vector<std::vector<std::string>> metric_labels_;

  // data representation
  KeyMap metric_map_;
  InvKeyMap inv_metric_map_;

  ResultsIndex index_;
};

void print_table(ResultsTable& table, std::ostream& ofstream = std::cout);

}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_RESULTS_H
