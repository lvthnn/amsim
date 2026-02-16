#include <amsim/output/estimator.h>
#include <amsim/params.h>
#include <amsim/state.h>
#include <amsim/utils.h>

namespace amsim {

std::vector<std::string> label_matrix(
    const std::vector<std::string>& row_labels,
    const std::vector<std::string>& col_labels,
    const std::optional<std::string>& row_suffix = std::nullopt,
    const std::optional<std::string>& col_suffix = std::nullopt) {
  std::size_t n_rows = row_labels.size();
  std::size_t n_cols = col_labels.size();
  std::vector<std::string> labels(n_rows * n_cols);

  for (std::size_t row = 0; row < n_rows; ++row)
    for (std::size_t col = 0; col < n_cols; ++col)
      labels[(row * n_rows) + col] =
          (row_labels[row] + row_suffix.value_or("")) +
          "::" + (col_labels[col] + col_suffix.value_or(""));

  return labels;
}

class EstimatorHeritability : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorHeritability(const Params& params)
      : PopulationEstimatorStrategy(
            params.sim.out_dir,
            "pheno_h2",
            params.pheno.names,
            params.pheno.n_pheno) {
    n_pheno_ = params.pheno.n_pheno;
  }

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno, 0) =
          state.pheno().comp_var(pheno, phenome::Component::Genetic) /
          state.pheno().comp_var(pheno, phenome::Component::Total);
  }

 private:
  std::size_t n_pheno_;
};

class EstimatorComponentMean : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentMean(
      const Params& params, phenome::Component type)
      : PopulationEstimatorStrategy(
            params.sim.out_dir,
            "pheno_" + phenome::to_string(type) + "_mean",
            params.pheno.names,
            params.pheno.n_pheno),
        type_(type),
        n_pheno_(params.pheno.n_pheno) {}

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno) = state.pheno().comp_mean(pheno, type_);
  }

 private:
  phenome::Component type_;
  std::size_t n_pheno_;
};

class EstimatorComponentVar : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentVar(
      const Params& params, phenome::Component type)
      : PopulationEstimatorStrategy(
            params.sim.out_dir,
            "pheno_" + phenome::to_string(type) + "_var",
            params.pheno.names,
            params.pheno.n_pheno),
        type_(type),
        n_pheno_(params.pheno.n_pheno) {}

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno) = state.pheno().comp_var(pheno, type_);
  }

 private:
  phenome::Component type_;
  std::size_t n_pheno_;
};

class EstimatorComponentCor : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorComponentCor(
      const Params& params,
      phenome::Component type_l,
      std::optional<phenome::Component> type_r)
      : PopulationEstimatorStrategy(
            params.sim.out_dir,
            "pheno_" + phenome::to_string(type_l) + "_" +
                phenome::to_string(type_r.value_or(type_l)) + "_cor",
            label_matrix(
                params.pheno.names,
                params.pheno.names,
                "_" + phenome::to_string(type_l),
                "_" + phenome::to_string(type_r.value_or(type_l))),
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
  phenome::Component type_l_;
  phenome::Component type_r_;
  Eigen::MatrixXd std_l_;
  Eigen::MatrixXd std_r_;
};

class EstimatorMateCor : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorMateCor(
      const Params& params, phenome::Component type)
      : PopulationEstimatorStrategy(
            params.sim.out_dir,
            "mate_" + phenome::to_string(type) + "_cor",
            label_matrix(
                params.pheno.names, params.pheno.names, "_male", "_female"),
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
  phenome::Component type_;
  std::size_t n_sex_;
  std::size_t n_pheno_;
  Eigen::MatrixXd std_male_;
  Eigen::MatrixXd std_female_;
};

PopulationEstimator PopulationHeritability() {
  return [](const Params& params) {
    return std::make_unique<EstimatorHeritability>(params);
  };
}

PopulationEstimator PopulationComponentMean(phenome::Component type) {
  return [type](const Params& params) {
    return std::make_unique<EstimatorComponentMean>(params, type);
  };
}

PopulationEstimator PopulationComponentVar(phenome::Component type) {
  return [type](const Params& params) {
    return std::make_unique<EstimatorComponentVar>(params, type);
  };
}

PopulationEstimator PopulationComponentCor(
    phenome::Component type_l,
    std::optional<phenome::Component> type_r) {
  return [type_l, type_r](const Params& params) {
    return std::make_unique<EstimatorComponentCor>(
        params, type_l, type_r);
  };
}

PopulationEstimator PopulationMateCor(phenome::Component type) {
  return [type](const Params& params) {
    return std::make_unique<EstimatorMateCor>(params, type);
  };
}

ComputePopulationEstimates::ComputePopulationEstimates(
    const Params& params,
    const PopulationEstimators& estimators,
    std::optional<std::size_t> rep_id)
    : rep_id_(std::move(rep_id)) {
  auto out_dir =
      (rep_id_.has_value())
          ? params.sim.out_dir / std::format("rep_{:03d}", rep_id_.value())
          : params.sim.out_dir;

  if (std::filesystem::exists(out_dir) && rep_id.has_value())
    throw std::runtime_error(
        "replicate directory " + out_dir.string() + " already exists!");
  std::filesystem::create_directory(out_dir);

  for (const auto& factory : estimators)
    estimators_.emplace_back(factory(params));
}

void ComputePopulationEstimates::operator()(const State& state) {
  for (auto& estimator : estimators_) (*estimator)(state);
}

}  // namespace amsim
