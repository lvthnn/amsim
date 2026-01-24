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
#include <array>

namespace amsim {

/// Number of summary statistics computed per metric
constexpr std::size_t NUM_STATS = 8;

/// @brief Map from string keys to indices
using KeyMap = std::unordered_map<std::string, std::size_t>;

/// @brief Map from indices to string keys
using InvKeyMap = std::unordered_map<std::size_t, std::string>;

/// @brief Table of summary statistics for a metric
///
/// ResultsTable stores summary statistics (mean, median, standard deviation,
/// standard error, and confidence intervals) across simulation replicates.
struct ResultsTable {
  ResultsTable() = default;

  /// @brief Add replicate values and compute summary statistics
  ///
  /// Computes mean, median, standard deviation, standard error, 95%
  /// confidence intervals, and 2.5%/97.5% quantiles.
  ///
  /// @param stream Vector of values from replicates
  void add(std::vector<double> stream) {
    std::size_t n_elem = stream.size();

    const double* ptr = stream.data();
    double mean = stats::mean(n_elem, ptr, 1);
    double sem = stats::sem(n_elem, ptr, 1);

    data[0].push_back(stats::mean(n_elem, ptr, 1));
    data[1].push_back(stats::quantile(0.5, n_elem, ptr, 1));
    data[2].push_back(stats::std(n_elem, ptr, 1, false));
    data[3].push_back(sem);
    data[4].push_back(mean - (1.96 * sem));
    data[5].push_back(mean + (1.96 * sem));
    data[6].push_back(stats::quantile(0.025, n_elem, ptr, 1));
    data[7].push_back(stats::quantile(0.975, n_elem, ptr, 1));
  }

  KeyMap label_map;                                 ///< Label to index map
  InvKeyMap inv_label_map;                          ///< Index to label map
  std::array<std::vector<double>, NUM_STATS> data;  ///< Summary statistics
};

/// @brief Vector of results tables indexed by metric
using ResultsIndex = std::vector<ResultsTable>;

/// @brief Manager for loading and summarising simulation results
///
/// SimulationResults loads metric data from simulation output files,
/// computes summary statistics across replicates, and saves results tables.
class SimulationResults {
 public:
  /// @brief Construct a SimulationResults instance
  ///
  /// @param out_dir Output directory containing simulation results
  /// @param n_replicates Optional number of replicates (inferred if not
  /// provided)
  /// @param metric_names Optional metric names (inferred if not provided)
  explicit SimulationResults(
      std::filesystem::path out_dir,
      std::optional<std::size_t> n_replicates = std::nullopt,
      std::optional<std::vector<std::string>> metric_names = std::nullopt);

  /// @brief Access results table for a specific metric
  ///
  /// @param metric Metric name
  /// @return ResultsTable for the metric
  ResultsTable operator()(const std::string& metric);

  /// @brief Compute summary statistics for all metrics
  void summarise();

  /// @brief Save results tables to disk
  ///
  /// @param metrics Optional subset of metrics to save
  /// @param out_dir Optional output directory override
  /// @param overwrite Whether to overwrite existing files
  void save(
      std::optional<std::vector<std::string>> metrics = std::nullopt,
      std::optional<std::filesystem::path> out_dir = std::nullopt,
      bool overwrite = false);

  /// @brief Return metric names
  /// @return Vector of metric names
  std::vector<std::string> metric_names() { return metric_names_; };

 private:
  /// @brief Infer number of replicates from output directory
  /// @param n_replicates Optional explicit replicate count
  void inferReplicates(std::optional<std::size_t> n_replicates = std::nullopt);

  /// @brief Infer metric names from output files
  /// @param metric_names Optional explicit metric names
  void inferMetrics(
      std::optional<std::vector<std::string>> metric_names = std::nullopt);

  /// @brief Encode string values as indices in maps
  /// @param values String values to encode
  /// @param map_out Output key-to-index map
  /// @param invmap_out Output index-to-key map
  static void encodeValues(
      const std::vector<std::string>& values,
      KeyMap& map_out,
      InvKeyMap& invmap_out);

  /// @brief Extract column labels from input streams
  /// @param labels Output vector of labels
  /// @param streams Input file streams
  static void getLabels(
      std::vector<std::string>& labels, std::vector<std::ifstream>& streams);

  /// @brief Compute summary statistics for a single metric
  /// @param metric Metric name
  void summariseMetric(const std::string& metric);

  std::filesystem::path out_dir_;  ///< Output directory path

  std::size_t n_replicates_ = 0;                 ///< Number of replicates
  std::vector<std::filesystem::path> rep_dirs_;  ///< Replicate directories

  std::size_t n_metrics_;                  ///< Number of metrics
  std::vector<std::string> metric_names_;  ///< Metric names
  std::vector<std::vector<std::string>>
      metric_labels_;  ///< Metric column labels

  KeyMap metric_map_;         ///< Metric name to index map
  InvKeyMap inv_metric_map_;  ///< Index to metric name map

  ResultsIndex index_;  ///< Results tables for all metrics
};

/// @brief Print a results table to an output stream
///
/// @param table ResultsTable to print
/// @param ofstream Output stream (default: stdout)
void print_table(ResultsTable& table, std::ostream& ofstream = std::cout);

}  // namespace amsim

#endif  // AMSIMCPP_SIMULATION_RESULTS_H
