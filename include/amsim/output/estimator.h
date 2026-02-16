#pragma once

#include <amsim/data/sample.h>
#include <amsim/params.h>
#include <amsim/state.h>
#include <amsim/utils.h>

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
PopulationEstimator PopulationComponentMean(
    phenome::Component type = phenome::Component::Total);
PopulationEstimator PopulationComponentVar(
    phenome::Component type = phenome::Component::Total);
PopulationEstimator PopulationComponentCor(
    phenome::Component type_l = phenome::Component::Total,
    std::optional<phenome::Component> type_r = std::nullopt);
PopulationEstimator PopulationMateCor(
    phenome::Component type = phenome::Component::Total);

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

// sample estimators

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
        const Params& params, const std::filesystem::path& sample_dir)>;

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
SampleEstimator<P> SampleMeanEstimator() {
  return [](const Params& params, const std::filesystem::path& sample_dir) {
    return std::make_unique<SampleMean<P>>(params, sample_dir);
  };
}

}  // namespace amsim
