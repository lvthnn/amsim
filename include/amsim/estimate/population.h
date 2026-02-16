#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>
#include <amsim/core/utils.h>

#include <amsim/sample/sampler.h>
#include <amsim/sample/proband.h>

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>

namespace amsim {

class PopulationEstimatorStrategy {
 public:
  PopulationEstimatorStrategy(
      const std::filesystem::path& out_dir,
      std::string name,
      std::vector<std::string> labels,
      std::size_t n_rows,
      std::size_t n_cols = 1)
      : name_(std::move(name)),
        labels_(std::move(labels)),
        n_rows_(n_rows),
        n_cols_(n_cols),
        stream_(out_dir / (name_ + ".tsv")),
        data_(n_rows_, n_cols_) {
    header();
  };

  virtual ~PopulationEstimatorStrategy() = default;
  virtual void compute(const State& state) = 0;

  void header() {
    stream_ << "gen\t";
    for (Eigen::Index el = 0; el < data_.size(); ++el) {
      stream_ << ((!labels_.empty()) ? labels_[el] : std::to_string(el));
      if (el < data_.size() - 1) stream_ << "\t";
    }
    stream_ << "\n";
  }

  void stream(std::size_t gen) {
    stream_ << std::to_string(gen);
    for (Eigen::Index el = 0; el < data_.size(); ++el)
      stream_ << "\t" + std::to_string(data_(el));
    stream_ << "\n";
  }

  void operator()(const State& state) {
    compute(state);
    stream(state.gen);
  }

 protected:
  std::string name_;
  std::vector<std::string> labels_;
  std::size_t n_rows_;
  std::size_t n_cols_;
  std::ofstream stream_;
  Eigen::MatrixXd data_;
};

using PopulationEstimator =
    std::function<std::unique_ptr<PopulationEstimatorStrategy>(
        const Params& params)>;

using PopulationEstimators = std::vector<PopulationEstimator>;

PopulationEstimator PopulationHeritability();
PopulationEstimator PopulationComponentMean(Component type = Component::Total);
PopulationEstimator PopulationComponentVar(Component type = Component::Total);
PopulationEstimator PopulationComponentCor(
    Component type_l = Component::Total,
    std::optional<Component> type_r = std::nullopt);
PopulationEstimator PopulationMateCor(Component type = Component::Total);

// macro class to manage all estimators simultaneously in simulation loop
class ComputePopulationEstimates {
 public:
  explicit ComputePopulationEstimates(
      const Params& params,
      const PopulationEstimators& estimators,
      std::optional<std::size_t> rep_id = std::nullopt);

  void operator()(const State& state);

 private:
  std::vector<std::unique_ptr<PopulationEstimatorStrategy>> estimators_;
  std::optional<std::size_t> rep_id_;
};

}  // namespace amsim
