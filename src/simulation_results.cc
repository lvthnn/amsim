#include <amsim/simulation_results.h>
#include <amsim/stats.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace amsim {

SimulationResults::SimulationResults(
    std::filesystem::path out_dir,
    std::optional<std::size_t> n_replicates,
    std::optional<std::vector<std::string>> metric_names)
    : out_dir_(std::move(out_dir)) {
  // ensure target directory exists
  if (!std::filesystem::exists(out_dir_))
    throw std::invalid_argument("specified output directory does not exist");

  // infer number of replicates and replicate directories if not suppied
  inferReplicates(n_replicates);
  inferMetrics(std::move(metric_names));
  encodeValues(metric_names_, metric_map_, inv_metric_map_);

  // resize class elements
  metric_labels_.resize(metric_names_.size());
  index_.resize(metric_names_.size());
}

ResultsTable SimulationResults::operator()(const std::string& metric) {
  auto it = metric_map_.find(metric);
  if (it == metric_map_.end())
    throw std::invalid_argument("metric " + metric + " not found");
  return index_[metric_map_[metric]];
}

void SimulationResults::encodeValues(
    const std::vector<std::string>& values,
    KeyMap& map_out,
    InvKeyMap& invmap_out) {
  map_out.reserve(values.size());

  std::size_t cnt = 0;
  for (const std::string& val : values) {
    map_out.insert({val, cnt});
    invmap_out.insert({cnt, val});
    ++cnt;
  }
}

void SimulationResults::inferReplicates(
    std::optional<std::size_t> n_replicates) {
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
}

void SimulationResults::inferMetrics(
    std::optional<std::vector<std::string>> metric_names) {
  if (!metric_names) {
    std::set<std::string> metrics_set;
    for (const auto& rep_dir : rep_dirs_) {
      std::set<std::string> rep_names;
      for (const auto& entry : std::filesystem::directory_iterator(rep_dir)) {
        if (entry.is_regular_file()) {
          std::string stem = entry.path().stem().string();
          rep_names.insert(stem);
        }
      }
      if (rep_dir == rep_dirs_.front()) {
        metrics_set = std::move(rep_names);
      } else {
        std::set<std::string> isect;
        std::ranges::set_intersection(
            metrics_set, rep_names, std::inserter(isect, isect.end()));
        metrics_set = std::move(isect);
      }
    }
    metric_names_ =
        std::vector<std::string>(metrics_set.begin(), metrics_set.end());
  } else {
    metric_names_ = *metric_names;
  }
  n_metrics_ = metric_names_.size();
}

void SimulationResults::getLabels(
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

void SimulationResults::summariseMetric(const std::string& metric) {
  std::vector<std::ifstream> streams;
  std::vector<double> stream_buf;
  std::vector<std::string> labels;

  streams.resize(n_replicates_);
  stream_buf.resize(n_replicates_);

  std::size_t metric_id = metric_map_[metric];

  // placerholder string token
  std::string dummy;
  std::string token;

  // create the ResultsTable to be used
  ResultsTable table;

  // initialise input stream vectors
  for (std::size_t rep = 0; rep < n_replicates_; ++rep) {
    std::filesystem::path file = rep_dirs_[rep] / (metric + ".tsv");
    streams[rep] = std::ifstream(file);
    if (!streams[rep].is_open())
      throw std::runtime_error(
          "could not open metric stream for file " + file.string());
  }

  // loop over each of the files
  getLabels(labels, streams);
  encodeValues(labels, table.label_map, table.inv_label_map);

  // summarise cells and push back into columns
  while (std::ranges::all_of(streams, [](std::ifstream& stream) {
    return stream.good() && stream.peek() != EOF;
  })) {
    // skip the first column (index) at the start of the row
    for (std::size_t rep = 0; rep < n_replicates_; ++rep)
      std::getline(streams[rep], dummy, '\t');

    // read and summarise columns
    for (std::size_t col = 0; col < labels.size(); ++col) {
      std::string label_name = table.inv_label_map[col];
      char delimiter = (col == labels.size() - 1) ? '\n' : '\t';

      for (std::size_t rep = 0; rep < n_replicates_; ++rep) {
        std::getline(streams[rep], token, delimiter);
        stream_buf[rep] = std::stod(token);
      }
      table.add(stream_buf);
    }
  }

  index_[metric_id] = std::move(table);
}

void SimulationResults::summarise() {
  for (const std::string& metric : metric_names_) summariseMetric(metric);
}

void SimulationResults::save(
    std::optional<std::vector<std::string>> metrics,
    std::optional<std::filesystem::path> out_dir,
    bool overwrite) {
  if (!metrics) metrics = metric_names_;
  if (!out_dir) out_dir = out_dir_;
  if (!std::filesystem::exists(*out_dir))
    std::filesystem::create_directory(*out_dir);

  for (const std::string& metric : *metrics) {
    std::filesystem::path metric_path = *out_dir / (metric + ".tsv");
    if (!overwrite && std::filesystem::exists(metric_path))
      throw std::invalid_argument(
          "file " + metric_path.string() + "already exists");
    std::ofstream metric_out(metric_path);

    if (!metric_out.is_open())
      throw std::invalid_argument(
          "could not open output stream for " + metric_path.string());

    print_table(index_[metric_map_[metric]], metric_out);
  }
}

void print_table(ResultsTable& table, std::ostream& ofstream) {
  std::size_t n_rows = table.data[0].size() / table.label_map.size();
  std::size_t n_cols = table.label_map.size();

  ofstream << "gen\tname\tmean\tmedian\tstddev\tstderr\tlower_ci95\tupper_"
              "ci95\tquant_025\tquant_975\n";

  ofstream << std::scientific;
  for (std::size_t col = 0; col < n_cols; ++col) {
    for (std::size_t row = 0; row < n_rows; ++row) {
      std::size_t idx = (row * n_cols) + col;
      ofstream << row + 1 << "\t" << table.inv_label_map[col] << "\t"
               << table.data[0][idx] << "\t" << table.data[1][idx] << "\t"
               << table.data[2][idx] << "\t" << table.data[3][idx] << "\t"
               << table.data[4][idx] << "\t" << table.data[5][idx] << "\t"
               << table.data[6][idx] << "\t" << table.data[7][idx] << "\n";
    }
  }
}

}  // namespace amsim
