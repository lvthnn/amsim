#pragma once

#include <amsim/core.h>
#include <amsim/sample/proband.h>

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>

namespace amsim {

template <Proband P>
class SampleEstimatorStrategy {
 public:
  SampleEstimatorStrategy(
      const std::filesystem::path& sample_dir,
      std::string name,
      std::vector<std::string> labels,
      std::size_t n_rows,
      std::size_t n_cols = 1)
      : name_(std::move(name)),
        labels_(std::move(labels)),
        n_rows_(n_rows),
        n_cols_(n_cols),
        stream_(sample_dir / (name_ + ".tsv")),
        data_(n_rows_, n_cols_) {
    header();
  };

  virtual ~SampleEstimatorStrategy() = default;
  virtual void compute(
      const Eigen::MatrixXd& phenotypes, const Eigen::MatrixXd& genotypes) = 0;

  std::string name() const { return name_; }

  void header() {
    stream_ << "gen";
    for (Eigen::Index el = 0; el < data_.size(); ++el)
      stream_ << "\t"
              << ((!labels_.empty()) ? labels_[el] : std::to_string(el));
    stream_ << "\n";
  }

  void stream(std::size_t gen) {
    stream_ << gen;
    for (Eigen::Index el = 0; el < data_.size(); ++el)
      stream_ << "\t" << data_(el);
    stream_ << "\n";
  }

  void operator()(
      std::size_t gen,
      const Eigen::MatrixXd& phenotypes,
      const Eigen::MatrixXd& genotypes) {
    compute(phenotypes, genotypes);
    stream(gen);
  }

 protected:
  std::string name_;
  std::vector<std::string> labels_;
  std::size_t n_rows_;
  std::size_t n_cols_;
  std::ofstream stream_;
  Eigen::MatrixXd data_;
};

template <Proband P>
using SampleEstimator =
    std::function<std::unique_ptr<SampleEstimatorStrategy<P>>(
        const Params& params,
        std::size_t n_probands,
        const std::filesystem::path& sample_dir)>;

template <Proband P>
using SampleEstimators = std::vector<SampleEstimator<P>>;

template <Proband P>
class SampleMean : public SampleEstimatorStrategy<P> {
 public:
  explicit SampleMean(
      const Params& params, const std::filesystem::path& sample_dir)
      : SampleEstimatorStrategy<P>(
            sample_dir,
            "sample_mean",
            params.pheno.names,
            params.pheno.n_pheno) {}

  void compute(
      const Eigen::MatrixXd& phenotypes,
      const Eigen::MatrixXd& /*genotypes*/) override {
    this->data_ = phenotypes.colwise().mean();
  }
};

template <Proband P>
class SampleVar : public SampleEstimatorStrategy<P> {
 public:
  explicit SampleVar(
      const Params& params, const std::filesystem::path& sample_dir)
      : SampleEstimatorStrategy<P>(
            sample_dir,
            "sample_var",
            params.pheno.names,
            params.pheno.n_pheno) {}

  void compute(
      const Eigen::MatrixXd& phenotypes,
      const Eigen::MatrixXd& /*genotypes*/) override {
    this->data_ = (phenotypes.rowwise() - phenotypes.colwise().mean())
                      .array()
                      .square()
                      .colwise()
                      .sum() /
                  (phenotypes.rows() - 1);
  }
};

template <Proband P>
class SampleMateCor : public SampleEstimatorStrategy<P> {
 public:
  explicit SampleMateCor(
      const Params& params,
      const std::filesystem::path& sample_dir,
      std::size_t n_probands)
      : SampleEstimatorStrategy<P>(
            sample_dir,
            "sample_mate_cor",
            utils::label_matrix(
                params.pheno.names, params.pheno.names, "_male", "_female"),
            params.pheno.n_pheno,
            params.pheno.n_pheno),
        n_probands_(n_probands),
        n_pheno_(params.pheno.n_pheno),
        std_male_(
            (n_probands_ * ProbandData<P>::ProbandSize) / 2,
            n_pheno_),
        std_female_(
            (n_probands_ * ProbandData<P>::ProbandSize) / 2,
            n_pheno_) {
    if constexpr (P == Proband::Individual) {
      throw std::runtime_error(
          "Unsupported proband type Proband::Individual for estimator "
          "SampleMateCor");
    }
  }

  // compute centres and standard deviations across columns
  void compute(
      const Eigen::MatrixXd& phenotypes,
      const Eigen::MatrixXd& /*genotypes*/) override {
    const std::size_t n_pairs = n_probands_ * 3;
    const std::size_t outer =
        n_probands_ * ProbandData<P>::ProbandSize;

    using MateMat = Eigen::Map<
        const Eigen::MatrixXd,
        Eigen::Unaligned,
        Eigen::Stride<Eigen::Dynamic, 2>>;

    MateMat males(
        phenotypes.data(),
        n_pairs,
        n_pheno_,
        {static_cast<Eigen::Index>(outer), 2});

    MateMat females(
        phenotypes.data() + 1,
        n_pairs,
        n_pheno_,
        {static_cast<Eigen::Index>(outer), 2});

    std_male_ = utils::standardise(males, false);
    std_female_ = utils::standardise(females, false);

    this->data_ = (1.0 / (n_pairs - 1)) * std_male_.transpose() * std_female_;
  }

 private:
  std::size_t n_probands_;
  std::size_t n_pheno_;
  Eigen::MatrixXd std_male_;
  Eigen::MatrixXd std_female_;
};

template <Proband P>
inline SampleEstimator<P> SampleMeanEstimator() {
  return [](const Params& params,
            std::size_t /*n_probands*/,
            const std::filesystem::path& sample_dir) {
    return std::make_unique<SampleMean<P>>(params, sample_dir);
  };
}

template <Proband P>
inline SampleEstimator<P> SampleVarEstimator() {
  return [](const Params& params,
            std::size_t /*n_probands*/,
            const std::filesystem::path& sample_dir) {
    return std::make_unique<SampleVar<P>>(params, sample_dir);
  };
}

template <Proband P>
inline SampleEstimator<P> SampleMateCorEstimator() {
  return [](const Params& params,
            std::size_t n_probands,
            const std::filesystem::path& sample_dir) {
    return std::make_unique<SampleMateCor<P>>(params, sample_dir, n_probands);
  };
}

}  // namespace amsim
