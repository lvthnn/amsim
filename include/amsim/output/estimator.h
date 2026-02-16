#pragma once

#include <amsim/params.h>
#include <amsim/state.h>
#include <amsim/utils.h>

#include <Eigen/Dense>
#include <fstream>

namespace amsim {

class EstimatorStrategy {
 public:
  EstimatorStrategy(
      std::string name,
      std::vector<std::string> labels,
      std::size_t n_rows,
      std::size_t n_cols = 1)
      : name_(std::move(name)),
        labels_(std::move(labels)),
        n_rows_(std::move(n_rows)),
        n_cols_(std::move(n_cols)),
        data_(n_rows_, n_cols_) {};

  virtual ~EstimatorStrategy() = default;
  virtual void compute(const State& state) = 0;

  std::string header() {
    std::string header = "gen\t";
    for (Eigen::Index el = 0; el < data_.size(); ++el) {
      header += (!labels_.empty()) ? labels_[el] : std::to_string(el);
      if (el < data_.size() - 1) header += "\t";
    }
    return header;
  }

  std::string stream(const State& state) {
    std::string res = std::to_string(state.gen) + "\t";
    for (Eigen::Index el = 0; el < data_.size(); ++el)
      res += "\t" + std::to_string(data_(el));

    return res;
  }

  std::string name() const { return name_; }

 protected:
  std::string name_;
  std::vector<std::string> labels_;
  std::size_t n_rows_;
  std::size_t n_cols_;
  Eigen::MatrixXd data_;
};

using Estimator = std::function<std::unique_ptr<EstimatorStrategy>(
    const Params& params)>;

using Estimators = std::vector<Estimator>;

Estimator Heritability();
Estimator ComponentMean(
    phenome::Component type = phenome::Component::Total);
Estimator ComponentVar(
    phenome::Component type = phenome::Component::Total);
Estimator ComponentCor(
    phenome::Component type_l = phenome::Component::Total,
    std::optional<phenome::Component> type_r = std::nullopt);
Estimator MateCor(phenome::Component type = phenome::Component::Total);

// macro class to manage all estimators simultaneously in simulation loop
class ComputeEstimates {
 public:
  explicit ComputeEstimates(
      const Params& params,
      const Estimators& estimators,
      std::optional<std::size_t> rep_id = std::nullopt);

  void headers();
  void operator()(const State& state);

 private:
  std::vector<std::unique_ptr<EstimatorStrategy>> estimators_;
  std::vector<std::ofstream> streams_;
  std::optional<std::size_t> rep_id_;
};

}  // namespace amsim
