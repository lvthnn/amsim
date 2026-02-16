#include <amsim/output/sample.h>
#include <amsim/state.h>

#include <numeric>

namespace amsim {

WeightFunction Uniform() {
  return [](const Eigen::MatrixXd& /*agg*/, Eigen::VectorXd& res) {
    res.setConstant(1.0);
  };
}

WeightFunction Logistic(const Eigen::VectorXd& effects) {
  return [effects](const Eigen::MatrixXd& agg, Eigen::VectorXd& res) {
    res = 1.0 / (1.0 + (-agg * effects).array().exp());
  };
}

// return the parity (parental / offspring) and index of the members of a
// proband that are subject to aggregation
template <Proband P>
void Sampler::Model<P>::fillAggregates(const State& state) {
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
void Sampler::Model<P>::extractProbands(const State& state) {
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
        std::size_t row = prob_id * n_members;
        std::size_t member = 0;

        if (of & Family::Son)
          phenotypes(row + member++, pheno) = pheno_vec(prob);
        if (of & Family::SonWife)
          phenotypes(row + member++, pheno) = pheno_vec(n_sex + matching[prob]);
        if (of & Family::Daughter)
          phenotypes(row + member++, pheno) =
              pheno_vec(n_sex + matching_par[prob]);
        if (of & Family::DaughterHusband)
          phenotypes(row + member++, pheno) =
              pheno_vec(inv_matching[matching_par[prob]]);
        if (of & Family::Father)
          phenotypes(row + member++, pheno) = pheno_par_vec(prob);
        if (of & Family::Mother)
          phenotypes(row + member++, pheno) =
              pheno_par_vec(n_sex + matching_par[prob]);
      }
    }

    else if constexpr (P == Proband::Mate) {
      for (std::size_t prob = 0; prob < n_probands; ++prob) {
        std::size_t prob_id = selected[prob];
        std::size_t row = prob_id * n_members;
        std::size_t member = 0;

        if (of & Mate::Husband)
          phenotypes(row + member++, pheno) = pheno_vec(prob);
        if (of & Mate::Wife)
          phenotypes(row + member++, pheno) = pheno_vec(n_sex + matching[prob]);
        if (of & Mate::HusbandFather)
          phenotypes(row + member++, pheno) = pheno_par_vec(prob);
        if (of & Mate::HusbandMother)
          phenotypes(row + member++, pheno) = pheno_par_vec(matching_par[prob]);
        if (of & Mate::WifeFather)
          phenotypes(row + member++, pheno) =
              pheno_par_vec(inv_matching[n_sex + matching[prob]]);
        if (of & Mate::WifeMother)
          phenotypes(row + member++, pheno) =
              pheno_par_vec(inv_matching_par[matching[prob]]);
      }
    }
  }
}

template <Proband P>
void Sampler::Model<P>::draw(const State& state) {
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
void Sampler::Model<P>::estimate(const State& state) {
  for (const auto& estimator : estimators)
    (*estimator)(state.gen, phenotypes, genotypes);
}

template <Proband P>
void Sampler::Model<P>::operator()(const State& state) {
  draw(state);
  estimate(state);
}

template struct Sampler::Model<Proband::Individual>;
template struct Sampler::Model<Proband::Mate>;
template struct Sampler::Model<Proband::Family>;

}  // namespace amsim
