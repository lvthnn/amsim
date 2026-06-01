#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>
#include <amsim/core/utils.h>
#include <amsim/io.h>
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
        data_(n_rows_, n_cols_) {
    Writer::create(data_, name_, row_labels_, col_labels_);
  }

  virtual ~PopulationEstimatorStrategy() = default;
  virtual void compute(const State& state) = 0;

  std::string name() const { return name_; }
  std::vector<std::string> row_labels() const { return row_labels_; }
  std::vector<std::string> col_labels() const { return col_labels_; }
  std::size_t n_rows() const { return n_rows_; }
  std::size_t n_cols() const { return n_cols_; }

  void operator()(const State& state) {
    compute(state);
    Writer::write(data_, name(), state.rep, state.gen);
  }

 protected:
  std::string name_;
  std::vector<std::string> row_labels_;
  std::vector<std::string> col_labels_;
  std::size_t n_rows_;
  std::size_t n_cols_;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> data_;
};

using PopulationEstimator =
    std::function<std::unique_ptr<PopulationEstimatorStrategy>(
        const Params& params)>;

using PopulationEstimators = std::vector<PopulationEstimator>;

namespace details {

class EstimatorGenotypeMean : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorGenotypeMean(const Params& params)
      : PopulationEstimatorStrategy("geno_mean", {}, {}, params.geno.n_loc) {};

  void compute(const State& state) override { data_ = state.geno().v_lmean(); }
};

class EstimatorGenotypeVar : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorGenotypeVar(const Params& params)
      : PopulationEstimatorStrategy("geno_var", {}, {}, params.geno.n_loc) {}

  void compute(const State& state) override { data_ = state.geno().v_lvar(); }
};

class EstimatorGenotypeMAF : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorGenotypeMAF(const Params& params)
      : PopulationEstimatorStrategy("geno_maf", {}, {}, params.geno.n_loc) {}

  void compute(const State& state) override { data_ = state.geno().v_lmaf(); }
};

class EstimatorGenotypeCov : public PopulationEstimatorStrategy {
 public:
  explicit EstimatorGenotypeCov(const Params& params)
      : PopulationEstimatorStrategy(
            "geno_cov", {}, {}, params.geno.n_loc, params.geno.n_loc) {
    n_loc_ = params.geno.n_loc;
  }

  void compute(const State& state) override {
    const GenoBuf& geno = state.geno();
    std::size_t n_words = geno.n_words();
    std::size_t n_ind = geno.n_ind();

    for (std::size_t loc1 = 0; loc1 < n_loc_; ++loc1) {
      // haplotypes of the first locus
      const uint64_t* h01 = geno.h0().rowptr(loc1);
      const uint64_t* h11 = geno.h1().rowptr(loc1);
      for (std::size_t loc2 = loc1; loc2 < n_loc_; ++loc2) {
        // haplotypes of the second locus
        const uint64_t* h02 = geno.h0().rowptr(loc2);
        const uint64_t* h12 = geno.h1().rowptr(loc2);
        std::size_t acc = 0;

        for (std::size_t word = 0; word < n_words; ++word)
          // can express the inner product as popcounts like so
          acc += __builtin_popcountll(h01[word] & h02[word]) +
                 __builtin_popcountll(h01[word] & h12[word]) +
                 __builtin_popcountll(h11[word] & h02[word]) +
                 __builtin_popcountll(h11[word] & h12[word]);

        data_(loc1, loc2) =
            (1.0 / n_ind) * acc - geno.v_lmean(loc1) * geno.v_lmean(loc2);
      }
    }

    data_.triangularView<Eigen::Lower>() =
        data_.transpose().triangularView<Eigen::Lower>();
  }

 private:
  std::size_t n_loc_;
};

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
            utils::vector_suffix(params.pheno.names, "_" + to_string(type_l)),
            utils::vector_suffix(
                params.pheno.names, "_" + to_string(type_r.value_or(type_l))),
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
            utils::vector_suffix(params.pheno.names, "_male"),
            utils::vector_suffix(params.pheno.names, "_female"),
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

inline PopulationEstimator PopulationGenotypeCov() {
  return [](const Params& params) {
    return std::make_unique<details::EstimatorGenotypeCov>(params);
  };
}

inline PopulationEstimator PopulationGenotypeMean() {
  return [](const Params& params) {
    return std::make_unique<details::EstimatorGenotypeMean>(params);
  };
}

inline PopulationEstimator PopulationGenotypeVar() {
  return [](const Params& params) {
    return std::make_unique<details::EstimatorGenotypeVar>(params);
  };
}

inline PopulationEstimator PopulationGenotypeMAF() {
  return [](const Params& params) {
    return std::make_unique<details::EstimatorGenotypeMAF>(params);
  };
}

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
