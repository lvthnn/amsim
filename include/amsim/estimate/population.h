#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>
#include <amsim/core/utils.h>
#include <amsim/sample/proband.h>

#include <Eigen/Dense>

namespace amsim {

class PopulationEstimatorStrategy {
 public:
  PopulationEstimatorStrategy(
      std::string name,
      std::vector<std::string> row_labels,
      std::vector<std::string> col_labels,
      std::size_t n_rows,
      std::size_t n_cols = 1)
      : name_(std::move(name)),
        row_labels_(std::move(row_labels)),
        col_labels_(std::move(col_labels)),
        n_rows_(n_rows),
        n_cols_(n_cols),
        data_(n_rows_, n_cols_) {}

  virtual ~PopulationEstimatorStrategy() = default;
  virtual void compute(const State& state) = 0;

  std::string name() const { return name_; }
  std::vector<std::string> row_labels() const { return row_labels_; }
  std::vector<std::string> col_labels() const { return col_labels_; }
  std::size_t n_rows() const { return n_rows_; }
  std::size_t n_cols() const { return n_cols_; }

  void operator()(const State& state) { compute(state); }

 protected:
  std::string name_;
  std::vector<std::string> row_labels_;
  std::vector<std::string> col_labels_;
  std::size_t n_rows_;
  std::size_t n_cols_;
  Eigen::MatrixXd data_;
};

using PopulationEstimator =
    std::function<std::unique_ptr<PopulationEstimatorStrategy>(
        const Params& params)>;

using PopulationEstimators = std::vector<PopulationEstimator>;

namespace details {

class EstimatorHeritability : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorHeritability(const Params& params)
      : PopulationEstimatorStrategy(
            "pheno_h2", {}, params.pheno.names, params.pheno.n_pheno) {
    n_pheno_ = params.pheno.n_pheno;
  }

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno) = state.pheno().comp_var(pheno, Component::Genetic) /
                     state.pheno().comp_var(pheno, Component::Total);
  }

 private:
  std::size_t n_pheno_;
};

class EstimatorComponentMean : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentMean(const Params& params, Component type)
      : PopulationEstimatorStrategy(
            "pheno_" + to_string(type) + "_mean",
            {},
            params.pheno.names,
            params.pheno.n_pheno),
        type_(type),
        n_pheno_(params.pheno.n_pheno) {}

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno) = state.pheno().comp_mean(pheno, type_);
  }

 private:
  Component type_;
  std::size_t n_pheno_;
};

class EstimatorComponentVar : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentVar(const Params& params, Component type)
      : PopulationEstimatorStrategy(
            "pheno_" + to_string(type) + "_var",
            {},
            params.pheno.names,
            params.pheno.n_pheno),
        type_(type),
        n_pheno_(params.pheno.n_pheno) {}

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno) = state.pheno().comp_var(pheno, type_);
  }

 private:
  Component type_;
  std::size_t n_pheno_;
};

class EstimatorComponentCor : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentCor(
      const Params& params, Component type_l, std::optional<Component> type_r)
      : PopulationEstimatorStrategy(
            "pheno_" + to_string(type_l) + "_" +
                to_string(type_r.value_or(type_l)) + "_cor",
            utils::label_vector(params.pheno.names, "_" + to_string(type_l)),
            utils::label_vector(
                params.pheno.names,
                "_" + to_string(type_r.value_or(type_l))),
            params.pheno.n_pheno,
            params.pheno.n_pheno),
        n_ind_(params.geno.n_ind),
        n_pheno_(params.pheno.n_pheno),
        type_l_(type_l),
        type_r_(type_r.value_or(type_l)),
        std_l_(n_ind_, n_pheno_),
        std_r_(n_ind_, n_pheno_) {}

  void compute(const State& state) override {
    std_l_ = utils::standardise(state.pheno()(type_l_));
    std_r_ = utils::standardise(state.pheno()(type_r_));
    data_ = (std_l_.transpose() * std_r_) / static_cast<double>(n_ind_);
  }

 private:
  std::size_t n_ind_;
  std::size_t n_pheno_;
  Component type_l_;
  Component type_r_;
  Eigen::MatrixXd std_l_;
  Eigen::MatrixXd std_r_;
};

class EstimatorMateCor : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorMateCor(const Params& params, Component type)
      : PopulationEstimatorStrategy(
            "mate_" + to_string(type) + "_cor",
            utils::label_vector(params.pheno.names, "_male"),
            utils::label_vector(params.pheno.names, "_female"),
            params.pheno.n_pheno,
            params.pheno.n_pheno),
        type_(type),
        n_sex_(params.geno.n_ind / 2),
        n_pheno_(params.pheno.n_pheno),
        std_male_(n_sex_, n_pheno_),
        std_female_(n_sex_, n_pheno_) {}

  void compute(const State& state) override {
    std_male_ = utils::standardise(state.pheno().male(type_));
    std_female_ = utils::standardise(state.pheno().female(type_));
    data_.setZero();

    for (std::size_t pair = 0; pair < n_sex_; ++pair)
      data_.noalias() += std_male_.row(pair).transpose() *
                         std_female_.row(state.matching()[pair]);
    data_ /= static_cast<double>(n_sex_);
  }

 private:
  Component type_;
  std::size_t n_sex_;
  std::size_t n_pheno_;
  Eigen::MatrixXd std_male_;
  Eigen::MatrixXd std_female_;
};

}  // namespace details

inline PopulationEstimator PopulationHeritability() {
  return [](const Params& params) {
    return std::make_unique<details::EstimatorHeritability>(params);
  };
}

inline PopulationEstimator PopulationComponentMean(
    Component type = Component::Total) {
  return [type](const Params& params) {
    return std::make_unique<details::EstimatorComponentMean>(params, type);
  };
}

inline PopulationEstimator PopulationComponentVar(
    Component type = Component::Total) {
  return [type](const Params& params) {
    return std::make_unique<details::EstimatorComponentVar>(params, type);
  };
}

inline PopulationEstimator PopulationComponentCor(
    Component type_l = Component::Total,
    std::optional<Component> type_r = std::nullopt) {
  return [type_l, type_r](const Params& params) {
    return std::make_unique<details::EstimatorComponentCor>(
        params, type_l, type_r);
  };
}

inline PopulationEstimator PopulationMateCor(
    Component type = Component::Total) {
  return [type](const Params& params) {
    return std::make_unique<details::EstimatorMateCor>(params, type);
  };
}

// macro class to manage all estimators simultaneously in simulation loop
class ComputePopulationEstimates {
 public:
  ComputePopulationEstimates(
      const Params& params, const PopulationEstimators& estimators) {
    for (const auto& factory : estimators)
      estimators_.emplace_back(factory(params));
  }

  void operator()(const State& state);

 private:
  std::vector<std::unique_ptr<PopulationEstimatorStrategy>> estimators_;
};

inline void ComputePopulationEstimates::operator()(const State& state) {
  for (auto& estimator : estimators_) (*estimator)(state);
}

}  // namespace amsim
