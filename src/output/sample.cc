#include <amsim/output/sample.h>
#include <amsim/state.h>

namespace amsim {

WeightFunction Uniform() {
  return [](const Eigen::VectorXd& /*agg*/, Eigen::VectorXd& res) {
    res.setConstant(1.0);
  };
}

WeightFunction Logistic(const Eigen::VectorXd& effects) {
  return [effects](const Eigen::VectorXd& agg, Eigen::VectorXd& res) {
    res = 1.0 / (1.0 + (-agg * effects).array().exp());
  };
}

// return the parity (parental / offspring) and index of the members of a
// proband that are subject to aggregation
template<Proband P>
std::array<std::size_t, 2> Sampler::Model<P>::getProbandMembers() {
  // for individuals:
  //   index of self = (parity, j)
  if (P == Proband::Individual) {}

  // for families:
  //   index of son = (parity, j)
  //   index of son's wife = (parity, matching[j])
  //   index of daughter = (parity, n_sex + matching[j]))
  //   index of daughter's husband = (parity, inv_matching[n_sex + matching[j]])
  //   index of dad = (1 - parity, j)
  //   index of mom = (1 - parity, n_sex + matching[j])
  if (P == Proband::Mate) {}

  // for mates:
  //   index of husband = (parity, j)
  //   index of husband's dad = (1 - parity, j)
  //   index of husband's mother = (1 - parity, n_sex + matching[j])
  //
  //   index of wife = (parity, matching[j])
  //   index of wife's dad = (1 - parity, inv_matching[n_sex + matching[j]])
  //   index of wife's mom = (1 - parity, matching[j])
  if (P == Proband::Mate) {}
}

template <Proband P>
void Sampler::Model<P>::draw(const State& state)  {
  // fill the aggregation buffer with phenotypes of the respective individuals
  // using get_proband_members()

  // run the aggregation buffer
}

template <Proband P>
void Sampler::Model<P>::estimate() {

}

}  // namespace amsim
