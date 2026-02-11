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

class EstimatorPhenotypeHeritability : public EstimatorImpl {
 public:
  explicit EstimatorPhenotypeHeritability(const Params& params)
      : EstimatorImpl("pheno_h2", params.pheno.names, params.pheno.n_pheno) {
    n_pheno_ = params.pheno.n_pheno;
  }

  void compute(const State& state) override {
    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
      data_(pheno, 0) =
          state.pheno().comp_var(pheno, phenome::ComponentType::GENETIC) /
          state.pheno().comp_var(pheno, phenome::ComponentType::TOTAL);
  }

 private:
  std::size_t n_pheno_;
};

class EstimatorPhenotypeComponentMean : public EstimatorImpl {
 public:
  explicit EstimatorPhenotypeComponentMean(
      const Params& params, phenome::ComponentType type)
      : EstimatorImpl(
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
  phenome::ComponentType type_;
  std::size_t n_pheno_;
};

class EstimatorPhenotypeComponentVar : public EstimatorImpl {
 public:
  explicit EstimatorPhenotypeComponentVar(
      const Params& params, phenome::ComponentType type)
      : EstimatorImpl(
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
  phenome::ComponentType type_;
  std::size_t n_pheno_;
};

class EstimatorPhenotypeComponentCor : public EstimatorImpl {
 public:
  explicit EstimatorPhenotypeComponentCor(
      const Params& params,
      phenome::ComponentType type_l,
      std::optional<phenome::ComponentType> type_r)
      : EstimatorImpl(
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
  phenome::ComponentType type_l_;
  phenome::ComponentType type_r_;
  Eigen::MatrixXd std_l_;
  Eigen::MatrixXd std_r_;
};

class EstimatorMateCorrelation : public EstimatorImpl {
 public:
  explicit EstimatorMateCorrelation(
      const Params& params, phenome::ComponentType type)
      : EstimatorImpl(
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
  phenome::ComponentType type_;
  std::size_t n_sex_;
  std::size_t n_pheno_;
  Eigen::MatrixXd std_male_;
  Eigen::MatrixXd std_female_;
};

Estimator PhenotypeHeritability() {
  return [](const Params& params) {
    return std::make_unique<EstimatorPhenotypeHeritability>(params);
  };
}

Estimator PhenotypeComponentMean(phenome::ComponentType type) {
  return [type](const Params& params) {
    return std::make_unique<EstimatorPhenotypeComponentMean>(params, type);
  };
}

Estimator PhenotypeComponentVar(phenome::ComponentType type) {
  return [type](const Params& params) {
    return std::make_unique<EstimatorPhenotypeComponentVar>(params, type);
  };
}

Estimator PhenotypeComponentCor(
    phenome::ComponentType type_l,
    std::optional<phenome::ComponentType> type_r) {
  return [type_l, type_r](const Params& params) {
    return std::make_unique<EstimatorPhenotypeComponentCor>(
        params, type_l, type_r);
  };
}

Estimator MateCorrelation(phenome::ComponentType type) {
  return [type](const Params& params) {
    return std::make_unique<EstimatorMateCorrelation>(params, type);
  };
}

ComputeEstimates::ComputeEstimates(
    const Params& params,
    const Estimators& estimators,
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

  std::size_t n_est = estimators.size();
  estimators_.resize(n_est);
  streams_.resize(n_est);

  for (std::size_t el = 0; el < n_est; ++el) {
    estimators_[el] = estimators[el](params);
    auto path = out_dir / (estimators_[el]->name() + ".tsv");

    streams_[el] = std::ofstream(path);

    if (!streams_[el].is_open())
      throw std::runtime_error(
          "estimator output " + path.string() + " could not be opened!");
  }
}

void ComputeEstimates::headers() {
  for (std::size_t el = 0; el < estimators_.size(); ++el)
    streams_[el] << estimators_[el]->header() << "\n";
}

void ComputeEstimates::operator()(const State& state) {
  for (std::size_t el = 0; el < estimators_.size(); ++el) {
    estimators_[el]->compute(state);
    streams_[el] << estimators_[el]->stream(state) << "\n";
  }
}

}  // namespace amsim
