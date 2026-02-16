#pragma once

#include <amsim/data/sample.h>
#include <amsim/output/estimator.h>
#include <amsim/params.h>
#include <amsim/state.h>

#include <filesystem>
#include <iostream>

namespace amsim {

WeightFunction Uniform();
WeightFunction Logistic(const Eigen::VectorXd& effects);

// user-facing specification struct
template <Proband P>
struct Sample {
  using ProbandEnum = ProbandData<P>::ProbandEnum;

  // selection parameters — which data, from whom, and what to do
  std::string name;
  std::size_t n_probands;
  std::optional<std::vector<std::string>> on;
  std::optional<std::vector<phenome::Component>> on_components;
  ProbandEnum of = ProbandData<P>::ProbandDefault;
  Aggregator agg = ProbandData<P>::AggDefault;
  bool decompress_genotypes = false;

  // weighting method — probability of selecting based on proband aggregate
  WeightFunction weighting;

  // estimators are declared here
  SampleEstimators<P> estimators;
};

using SampleSpec = std::variant<
    Sample<Proband::Individual>,
    Sample<Proband::Mate>,
    Sample<Proband::Family>>;

class Sampler {
 public:
  template <Proband P>
  explicit Sampler(Sample<P> sample, const Params& params)
      : self_(std::make_unique<Model<P>>(std::move(sample), params)){};

  void operator()(const State& state) {
    (*self_)(state);
  }

 private:
  // type-erased computing interface
  struct Concept {
    virtual ~Concept() = default;
    virtual void draw(const State& state) = 0;
    virtual void estimate(const State& state) = 0;
    virtual void operator()(const State& state) = 0;
  };

  // implementation strategy
  template <Proband P>
  struct Model : Concept {
    explicit Model(Sample<P> sample, const Params& params)
        : name(sample.name),
          n_probands(sample.n_probands),
          of(sample.of),
          agg(std::move(sample.agg)),
          weighting(std::move(sample.weighting)),
          decompress_genotypes(sample.decompress_genotypes),
          n_probands_total(
              P == Proband::Individual ? params.geno.n_ind
                                       : params.geno.n_ind / 2),
          n_sex(params.geno.n_ind / 2),
          n_members(__builtin_popcountll(static_cast<uint8_t>(sample.of))),
          n_pheno(params.pheno.n_pheno),
          n_loc(params.geno.n_loc),
          n_on(sample.on.has_value() ? sample.on.value().size() : 1),
          on_components(
              sample.on_components.has_value()
                  ? std::move(sample.on_components.value())
                  : std::vector<phenome::Component>(
                        n_on, phenome::Component::Total)),
          members(n_probands_total * n_members, n_on),
          aggregates(n_probands_total, n_on),
          unif(n_probands_total),
          keys(n_probands_total),
          phenotypes(sample.n_probands * ProbandData<P>::ProbandSize, n_pheno),
          selected(n_probands_total) {
      if (sample.on.has_value())
        for (const auto& pheno_name : sample.on.value())
          on_indices.push_back(params.pheno.pheno_ids.at(pheno_name));

      auto sample_dir = params.sim.out_dir / name;
      std::filesystem::create_directories(sample_dir);
      for (const auto& estimator : sample.estimators)
        estimators.emplace_back(estimator(params, sample_dir));
    };

    // fields obtained from sampler specification
    std::string name;
    std::size_t n_probands;
    ProbandData<P>::ProbandEnum of;
    Aggregator agg;
    WeightFunction weighting;
    std::vector<std::unique_ptr<SampleEstimatorStrategy<P>>> estimators;
    bool decompress_genotypes;

    // derived dimensions
    std::size_t n_probands_total;
    std::size_t n_sex;
    std::size_t n_members;
    std::size_t n_pheno;
    std::size_t n_loc;
    std::size_t n_on;

    // phenotype indices and components for sampling
    std::vector<std::size_t> on_indices;
    std::vector<phenome::Component> on_components;

    // buffers to store sampling keys and data
    Eigen::MatrixXd members;
    Eigen::MatrixXd aggregates;
    Eigen::VectorXd unif;
    Eigen::VectorXd keys;
    Eigen::MatrixXd phenotypes;
    Eigen::MatrixXd genotypes;
    std::vector<std::size_t> selected;

    void fillAggregates(const State& state);
    void extractProbands(const State& state);

    void draw(const State& state) override;
    void estimate(const State& state) override;
    void operator()(const State& state) override;
  };

  std::unique_ptr<Concept> self_;
};

}  // namespace amsim

// The above sampler allows:
//
// Sample<Proband::Individual> random_sample{
//   .n_probands = 12000,
//   .weighting = Uniform(),
//   .estimators = {Heritability(), ComponentMean()}
// };
//
// Sample<Proband::Family> participation_bias{
//   .n_probands = 150,
//   .on = {"educational_attainment"},
//   .on_components = {Component::Genetic},
//   .of = Family::Parents,
//   .weighting = Logistic(Eigen::VectorXd{2.0}),
//   .estimators = {Heritability(), GWAS()}
// };
//
// Sample<Proband::Mate> am_pb{
//   .n_probands = 150,
//   .on = {"educational_attainment"},
//   .on_components = {Component::Total},
//   .of = Mate::Husband,
//   .weighting = Logistic(0.5),
//   .estimators = {MateCorrelation()}
// };
