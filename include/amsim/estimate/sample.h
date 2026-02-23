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
  }

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
        std_male_((n_probands_ * ProbandData<P>::ProbandSize) / 2, n_pheno_),
        std_female_((n_probands_ * ProbandData<P>::ProbandSize) / 2, n_pheno_) {
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
    const std::size_t outer = n_probands_ * ProbandData<P>::ProbandSize;

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
class SampleGWAS : public SampleEstimatorStrategy<P> {
 public:
  explicit SampleGWAS(
      const Params& params, const std::filesystem::path& sample_dir)
      : SampleEstimatorStrategy<P>(params, sample_dir) {
    utils::check_plink2();
  }

 private:
};

template <Proband P>
class SampleHasemanElston : public SampleEstimatorStrategy<P> {
 public:
  explicit SampleHasemanElston(
      const Params& params, const std::filesystem::path& sample_dir)
      : SampleHasemanElston<P>(params, sample_dir) {
    utils::check_gcta64();
  }

 private:
};

template <Proband P>
class SampleGREML : public SampleEstimatorStrategy<P> {
 public:
  explicit SampleGREML(
      const Params& params, const std::filesystem::path& sample_dir)
      : SampleGREML<P>(params, sample_dir) {
    utils::check_gcta64();
  }

 private:
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
inline SampleEstimator<P> SampleCovEstimator() {}

template <Proband P>
inline SampleEstimator<P> SampleMateCorEstimator() {
  return [](const Params& params,
            std::size_t n_probands,
            const std::filesystem::path& sample_dir) {
    return std::make_unique<SampleMateCor<P>>(params, sample_dir, n_probands);
  };
}

// invoke PLINK2
// should support PCA correction for ancestry
// - both a boolean flag whether to do that,
// - then flags which allow controlling how many PCA components
template <Proband P>
inline SampleEstimator<P> SampleGWASEstimator() {
  // whoops, this is actually meant for the compute function
  //
  // std::string gcta_string = std::format(...);
  // std::system(gcta_string);
  //
  // verify call worked
  //
  // extract values from file, and compute:
  //   - true positive rate
  //   - false positive rate
  //   - bias / l2 / linfty norms
  //   - more?
  // write to HDF5
  //
  // delete the intermediate file
  return [](const Params& params,
            std::size_t n_probands,
            const std::filesystem::path& sample_dir) {
    return std::make_unique<SampleGWAS<P>>(params, n_probands, sample_dir);
  };
}

// invoke GCTA
template <Proband P>
inline SampleEstimator<P> SampleHasemanElstonEstimator() {
  // whoops, this is actually meant for the compute function
  //
  // std::string gcta_string = std::format(...);
  // std::system(gcta_string);
  //
  // verify call worked
  //
  // extract the values from the file and write it to HDF5
  //
  // delete the intermediate file
}

// invoke GCTA
template <Proband P>
inline SampleEstimator<P> SampleGREMLEstimator() {
  // whoops, this is actually meant for the compute function
  //
  // support pca correction? likely unnecessary since there is no ancestry
  //
  // std::string gcta_string = std::format(...);
  // std::system(gcta_string);
  //
  // verify call worked
  //
  // extract the GREML estimate and write it to HDF5
  //
  // delete the intermediate files
}

// invoke PLINK and compute PCA
template <Proband P>
inline SampleEstimator<P> SamplePCAEstimator() {
  // is this unnecessary?
  //
  // would be very interesting to see how PCA vectors and eigenvalues
  // change with time
}

}  // namespace amsim
