#pragma once

#include <amsim/core.h>
#include <amsim/estimate.h>
#include <amsim/sample.h>

#include <filesystem>

namespace amsim {

inline WeightFunction Uniform() {
  return [](const Eigen::MatrixXd& /*agg*/, Eigen::VectorXd& res) {
    res.setConstant(1.0);
  };
}

inline WeightFunction Logistic(const Eigen::VectorXd& effects) {
  return [effects](const Eigen::MatrixXd& agg, Eigen::VectorXd& res) {
    res = 1.0 / (1.0 + (-agg * effects).array().exp());
  };
}

// user-facing specification struct
template <Proband P>
struct Sample {
  using ProbandEnum = ProbandData<P>::ProbandEnum;

  // selection parameters — which data, from whom, and what to do
  std::string name;
  std::size_t n_probands;
  std::optional<std::vector<std::string>> on;
  std::optional<std::vector<Component>> on_components;
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

  void operator()(const State& state) { (*self_)(state); }

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
                  : std::vector<Component>(n_on, Component::Total)),
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
    std::vector<Component> on_components;

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

template <Proband P>
inline void Sampler::Model<P>::fillAggregates(const State& state) {
  if (on_indices.empty()) return;
  const auto& matching = state.matching();
  const auto& matching_par = state.matching_par();
  const auto& inv_matching = state.inv_matching();
  const auto& inv_matching_par = state.inv_matching_par();

  for (std::size_t on = 0; on < n_on; ++on) {
    const auto& pheno_on = state.pheno()(on_indices[on], on_components[on]);
    const auto& pheno_on_par =
        state.pheno_par()(on_indices[on], on_components[on]);

    if constexpr (P == Proband::Individual)
      members.col(on) = pheno_on;

    else if constexpr (P == Proband::Family) {
      for (std::size_t prob = 0; prob < n_probands_total; ++prob) {
        std::size_t member = 0;
        std::size_t row = prob * n_members;

        if (of & Family::Son) members(row + member++, on) = pheno_on(prob);

        if (of & Family::SonWife)
          members(row + member++, on) = pheno_on(n_sex + matching[prob]);

        if (of & Family::Daughter)
          members(row + member++, on) = pheno_on(n_sex + matching_par[prob]);

        if (of & Family::DaughterHusband)
          members(row + member++, on) =
              pheno_on(inv_matching[matching_par[prob]]);

        if (of & Family::Father)
          members(row + member++, on) = pheno_on_par(prob);

        if (of & Family::Mother)
          members(row + member++, on) =
              pheno_on_par(n_sex + matching_par[prob]);
      }
    }

    else if constexpr (P == Proband::Mate) {
      for (std::size_t prob = 0; prob < n_probands_total; ++prob) {
        std::size_t member = 0;
        std::size_t row = prob * n_members;

        if (of & Mate::Husband) members(row + member++, on) = pheno_on(prob);

        if (of & Mate::Wife)
          members(row + member++, on) = pheno_on(n_sex + matching[prob]);

        if (of & Mate::HusbandFather)
          members(row + member++, on) = pheno_on_par(prob);

        if (of & Mate::HusbandMother)
          members(row + member++, on) =
              pheno_on_par(n_sex + matching_par[prob]);

        if (of & Mate::WifeFather)
          members(row + member++, on) =
              pheno_on_par(inv_matching[n_sex + matching[prob]]);

        if (of & Mate::WifeMother)
          members(row + member++, on) =
              pheno_on_par(inv_matching_par[matching[prob]]);
      }
    }
  }

  // aggregate each of the probands
  for (std::size_t prob = 0; prob < n_probands_total; ++prob) {
    auto proband = members.middleRows(prob * n_members, n_members);

    if (agg == Aggregator::Mean)
      aggregates.row(prob) = proband.colwise().mean();
    else if (agg == Aggregator::Max)
      aggregates.row(prob) = proband.colwise().maxCoeff();
    else if (agg == Aggregator::Min)
      aggregates.row(prob) = proband.colwise().minCoeff();
    else if (agg == Aggregator::Identity)
      aggregates.row(prob) = proband.row(0);
  }
}

template <Proband P>
inline void Sampler::Model<P>::extractProbands(const State& state) {
  const auto& matching = state.matching();
  const auto& matching_par = state.matching_par();
  const auto& inv_matching = state.inv_matching();
  const auto& inv_matching_par = state.inv_matching_par();

  for (std::size_t pheno = 0; pheno < n_pheno; ++pheno) {
    const auto& pheno_vec = state.pheno()(pheno);
    const auto& pheno_par_vec = state.pheno_par()(pheno);

    if constexpr (P == Proband::Individual)
      for (std::size_t prob = 0; prob < n_probands; ++prob)
        phenotypes(prob, pheno) = pheno_vec(selected[prob]);

    else if constexpr (P == Proband::Family) {
      for (std::size_t prob = 0; prob < n_probands; ++prob) {
        std::size_t prob_id = selected[prob];
        std::size_t row = prob_id * ProbandData<P>::ProbandSize;

        phenotypes(row, pheno) = pheno_vec(prob_id);
        phenotypes(row + 1, pheno) = pheno_vec(n_sex + matching[prob_id]);
        phenotypes(row + 2, pheno) = pheno_vec(n_sex + matching_par[prob_id]);
        phenotypes(row + 3, pheno) =
            pheno_vec(inv_matching[matching_par[prob_id]]);
        phenotypes(row + 4, pheno) = pheno_par_vec(prob_id);
        phenotypes(row + 5, pheno) =
            pheno_par_vec(n_sex + matching_par[prob_id]);
      }
    }

    else if constexpr (P == Proband::Mate) {
      for (std::size_t prob = 0; prob < n_probands; ++prob) {
        std::size_t prob_id = selected[prob];
        std::size_t row = prob_id * ProbandData<P>::ProbandSize;

        phenotypes(row, pheno) = pheno_vec(prob_id);
        phenotypes(row + 1, pheno) = pheno_vec(n_sex + matching[prob_id]);
        phenotypes(row + 2, pheno) = pheno_par_vec(prob_id);
        phenotypes(row + 3, pheno) = pheno_par_vec(matching_par[prob_id]);
        phenotypes(row + 4, pheno) =
            pheno_par_vec(inv_matching[n_sex + matching[prob_id]]);
        phenotypes(row + 5, pheno) =
            pheno_par_vec(inv_matching_par[matching[prob_id]]);
      }
    }
  }
}

template <Proband P>
inline void Sampler::Model<P>::draw(const State& state) {
  // fill the aggregation buffer with phenotypes of the respective members
  fillAggregates(state);

  // generate uniform weights for Efraimidis-Spirakis
  rng::UniformRange::fill(unif.data(), n_probands_total);

  // run the aggregation buffer
  weighting(aggregates, keys);
  keys = unif.array().pow(1.0 / keys.array());

  std::iota(selected.begin(), selected.end(), 0);
  std::nth_element(
      selected.begin(),
      selected.begin() + n_probands,
      selected.end(),
      [&](auto a, auto b) { return keys(a) > keys(b); });

  // extract selected proband data
  extractProbands(state);
}

template <Proband P>
inline void Sampler::Model<P>::estimate(const State& state) {
  for (const auto& estimator : estimators)
    (*estimator)(state.gen, phenotypes, genotypes);
}

template <Proband P>
inline void Sampler::Model<P>::operator()(const State& state) {
  draw(state);
  estimate(state);
}

template struct Sampler::Model<Proband::Individual>;
template struct Sampler::Model<Proband::Mate>;
template struct Sampler::Model<Proband::Family>;

}  // namespace amsim
