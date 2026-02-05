#pragma once

#include <amsim/params.h>
#include <amsim/state.h>

#include <Eigen/Dense>
#include <fstream>

#include <amsim/utils.h>

namespace amsim {

class EstimatorImpl {
 public:
  EstimatorImpl(
      std::string name,
      std::vector<std::string> labels,
      std::size_t n_rows,
      std::size_t n_cols = 1)
      : name_(std::move(name)),
        labels_(std::move(labels)),
        n_rows_(std::move(n_rows)),
        n_cols_(std::move(n_cols)),
        data_(n_rows_, n_cols_) {};

  virtual ~EstimatorImpl() = default;
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
    for (std::size_t row = 0; row < n_rows_; ++row)
      for (std::size_t col = 0; col < n_cols_; ++col)
        res += std::to_string(data_(row, col)) +
               (row * col < ((n_rows_ * n_cols_) - 1) ? "\t" : "");

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

using Estimator =
    std::function<std::unique_ptr<EstimatorImpl>(const Params& params)>;

using Estimators = std::vector<Estimator>;

Estimator PhenoHeritability();
Estimator PhenoComponentMean(
    phenome::ComponentType type = phenome::ComponentType::TOTAL);
Estimator PhenoComponentVar(
    phenome::ComponentType type = phenome::ComponentType::TOTAL);
Estimator PhenoComponentCor(
    phenome::ComponentType type_l = phenome::ComponentType::TOTAL,
    std::optional<phenome::ComponentType> type_r = std::nullopt);
Estimator MateCorrelation(
    phenome::ComponentType type = phenome::ComponentType::TOTAL);

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
  std::vector<std::unique_ptr<EstimatorImpl>> estimators_;
  std::vector<std::ofstream> streams_;
  std::optional<std::size_t> rep_id_;
};

}  // namespace amsim
