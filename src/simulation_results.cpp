#include <amsim/simulation_results.h>
#include <amsim/stats.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace amsim {
CellStats::CellStats(std::vector<double> stream)
    : n_elem(stream.size()),
      mean(stats::mean(n_elem, stream.data(), 1)),
      median(stats::quantile(0.5, n_elem, stream.data(), 1)),
      std(std::sqrt(stats::var(n_elem, stream.data(), 1, false))),
      sem(std / std::sqrt(n_elem)),
      lci(mean - 1.96 * sem),
      uci(mean + 1.96 * sem),
      q025(stats::quantile(0.025, n_elem, stream.data(), 1)),
      q975(stats::quantile(0.975, n_elem, stream.data(), 1)) {}

SimulationResults::SimulationResults(
    std::filesystem::path out_dir,
    std::optional<std::size_t> n_replicates,
    std::optional<std::set<std::string>> metric_names)
    : out_dir_(std::move(out_dir)) {
  if (!std::filesystem::exists(out_dir_)) {
    throw std::invalid_argument("specified output directory does not exist");
  }

  // infer number of replicates and replicate directories if not suppied
  if (!n_replicates) {
    for (const auto& dir_entry :
         std::filesystem::directory_iterator(out_dir_)) {
      if (dir_entry.is_directory()) {
        const std::string stem = dir_entry.path().stem().string();
        const std::regex rep_regex("rep_\\d{3}");
        if (std::regex_match(stem, rep_regex)) {
          ++n_replicates_;
          rep_dirs_.push_back(dir_entry.path());
        }
      }
    }
  } else {
    n_replicates_ = *n_replicates;
  }

  // infer metric names if not supplied
  if (!metric_names) {
    std::set<std::string> metric_names;

    for (const auto& rep_dir : rep_dirs_) {
      std::set<std::string> rep_names;
      for (const auto& entry : std::filesystem::directory_iterator(rep_dir)) {
        if (entry.is_regular_file()) {
          rep_names.insert(entry.path().stem().string());
        }
      }
      if (rep_dir == rep_dirs_.front()) {
        metric_names = rep_names;
      } else {
        std::set_intersection(
            metric_names.begin(),
            metric_names.end(),
            rep_names.begin(),
            rep_names.end(),
            std::inserter(metric_names_, metric_names_.end()));
      }
    }
  } else {
    metric_names_ = *metric_names;
  }

  // set up encoding of metric names as integers
  std::size_t map_names = 0;
  for (const auto& metric : metric_names_) {
    metric_map_.insert({metric, map_names});
    ++map_names;
  }

  // resize the outermost indexer
  table_.resize(metric_names_.size());
  label_map_.resize(metric_names_.size());
  label_names_.resize(metric_names_.size());
}

void get_labels(
    std::vector<std::string>& labels, std::vector<std::ifstream>& streams) {
  std::string header_line;
  std::string token;
  std::string dummy;

  // read the header, get labels, and skip index
  std::getline(streams[0], header_line);
  std::stringstream ss(header_line);
  std::getline(ss, dummy, '\t');

  // push labels back to the vector
  while (std::getline(ss, token, '\t')) {
    if (!token.empty() && token.back() == '\n') token.pop_back();
    labels.push_back(token);
  }

  // throw away header in other streams
  for (std::size_t rep = 1; rep < streams.size(); ++rep) {
    std::getline(streams[rep], dummy);
  }
}

void SimulationResults::summarise_metric_(std::string metric_name) {
  std::vector<std::ifstream> streams(n_replicates_);
  std::vector<double> stream_buf(n_replicates_);
  std::vector<std::string> labels;

  std::string dummy, token;

  // initialise each replicate stream
  for (std::size_t rep = 0; rep < n_replicates_; ++rep) {
    std::filesystem::path file = rep_dirs_[rep] / (metric_name + ".tsv");
    streams[rep] = std::ifstream(file);
  }

  get_labels(labels, streams);

  label_names_[metric_map_[metric_name]].resize(labels.size());
  label_names_[metric_map_[metric_name]] = labels;

  KeyMap& metric_label_map = label_map_[metric_map_[metric_name]];
  metric_label_map.reserve(labels.size());
  for (std::size_t it = 0; it < labels.size(); it++) {
    metric_label_map.insert({labels[it], it});
  }

  table_[metric_map_[metric_name]].resize(labels.size());

  std::size_t row = 0;
  while (std::all_of(streams.begin(), streams.end(), [](std::ifstream& stream) {
    return stream.good() && stream.peek() != EOF;
  })) {
    // skip the first column (index) at the start of the row
    for (std::size_t rep = 0; rep < n_replicates_; ++rep) {
      std::getline(streams[rep], dummy, '\t');
    }

    // read and summarise columns
    for (std::size_t col = 0; col < labels.size(); ++col) {
      std::string label_name = labels[col];
      char delimiter = (col == labels.size() - 1) ? '\n' : '\t';

      for (std::size_t rep = 0; rep < n_replicates_; ++rep) {
        std::getline(streams[rep], token, delimiter);
        stream_buf[rep] = std::stod(token);
      }

      CellStats cell_stats(stream_buf);
      std::size_t metric_id = metric_map_[metric_name];
      std::size_t label_id = metric_label_map[label_name];

      table_[metric_id][label_id].push_back(cell_stats);
    }
    ++row;
  }
}

void SimulationResults::summarise() {
  for (const auto& metric : metric_names_) {
    summarise_metric_(metric);
  }
}

void SimulationResults::print(std::string metric_name) {
  std::vector<std::vector<CellStats>> metric_data = table_[metric_map_[metric_name]];
  std::cout << "gen\tname\tmean\tmedian\tstd\tsem\tq025\tq975\tlci\tuci\n";

  for (std::size_t row = 0; row < metric_data.size(); row++) {
    std::vector<CellStats> stats = metric_data[row];
    for (std::size_t col = 0; col < stats.size(); col++) {
      CellStats cell_stat = stats[col];
      std::cout << col + 1 << "\t"
                << label_names_[metric_map_[metric_name]][row] << "\t"
                << cell_stat.mean << "\t"
                << cell_stat.median << "\t"
                << cell_stat.std << "\t"
                << cell_stat.sem << "\t"
                << cell_stat.q025 << "\t"
                << cell_stat.q975 << "\t"
                << cell_stat.lci << "\t"
                << cell_stat.uci << "\n";
    }
  }
}

}  // namespace amsim
