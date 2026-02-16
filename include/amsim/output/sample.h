#pragma once

#include <amsim/output/estimator.h>
#include <amsim/data/sample.h>
#include <amsim/params.h>
#include <amsim/state.h>

namespace amsim {

WeightFunction Uniform();
WeightFunction Logistic(Eigen::VectorXd effects);

// user-facing specification struct
template <Proband P>
struct Sample {
  using ProbandEnum = ProbandData<P>::ProbandEnum;

  // selection parameters — which data, from whom, and what to do
  std::size_t n_probands;
  std::optional<std::vector<std::string>> on;
  std::optional<std::vector<phenome::Component>> on_components;
  ProbandEnum of = ProbandData<P>::ProbandDefault;
  Aggregator agg = ProbandData<P>::AggDefault;
  bool decompress_genotypes = false;

  // weighting method — probability of selecting based on proband aggregate
  WeightFunction weighting;
  
  // estimators are declared here
  Estimators estimators;
};

class Sampler {
 public:
  template <Proband P>
  explicit Sampler(Sample<P> sample, const Params& params)
    : n_ind_(params.geno.n_ind),
      n_loc_(params.geno.n_loc),
      n_pheno_(params.pheno.n_pheno),
      self_(std::make_unique<Model<P>>(std::move(sample))) {};

 private:
  struct Concept {
    virtual ~Concept() = default;
    virtual void draw(const State& state) = 0; 
    virtual void estimate() = 0;
  };

  template <Proband P>
  struct Model : Concept {
    explicit Model(Sample<P> sample_) : sample(std::move(sample_)) {};

    Sample<P> sample;
    std::vector<std::size_t> idx; 
    Eigen::VectorXd keys;
    Eigen::MatrixXd phenotypes;
    std::array<std::size_t, 2> getProbandMembers();
    void draw(const State& state) override;
    void estimate() override;
  };

  std::size_t n_ind_;
  std::size_t n_loc_;
  std::size_t n_pheno_;
  std::unique_ptr<Concept> self_;
};

template struct Sampler::Model<Proband::Individual>;
template struct Sampler::Model<Proband::Family>;
template struct Sampler::Model<Proband::Mate>;

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
