#ifndef AMSIMCPP_MATING_H
#define AMSIMCPP_MATING_H

#include <amsim/mating_type.h>
#include <amsim/phenotype.h>
#include <amsim/rng.h>

#include <cstddef>
#include <random>
#include <vector>

namespace amsim {
class MatingModel {
 public:
  explicit MatingModel(MatingType type, std::size_t n_sex)
      : type_(type), g_(std::random_device{}()), n_sex_(n_sex) {}

  virtual ~MatingModel() = default;
  virtual std::vector<std::size_t> match() = 0;

 protected:
  MatingType type_;
  std::mt19937 g_;
  std::size_t n_sex_;

  std::vector<std::size_t> randState();
};

class AssortativeModel : public MatingModel {
 public:
  AssortativeModel(
      const PhenotypeList& phenotypes,
      std::vector<double> cor,
      std::size_t n_sex,
      const rng::Xoshiro256ss& rng,
      std::size_t n_itr = 2e6,
      double temp_init = 0.50,
      double temp_decay = 0.9999,
      double tol_inf = 1e-7);

  void display_cor();
  void init_state();
  std::vector<std::size_t> match() override;
  void update(const PhenotypeList& phenotypes);

 private:
  std::vector<const double*> ptr_tot_;

  const std::vector<double> cor_;
  const std::size_t n_pheno_;
  const std::size_t n_sex_;
  const double tol_inf_;
  const std::size_t n_itr_;
  const double temp_init_;
  const double temp_decay_;

  std::vector<double> male_;
  std::vector<double> female_;

  void arrange();

  std::vector<double> computeCor();
  std::vector<double> computeDelta(std::size_t i0, std::size_t i1);

  double computeDiffEnergy(
      const std::vector<double>& cur,
      const std::vector<double>& target,
      const std::vector<double>& delta) const;

  rng::NormalPolar fuzz_;
  rng::UniformIntRange swap_;
  rng::UniformRange acc_;

 public:
  std::vector<double> cor_S;
  std::vector<double> cor_U;
  std::vector<double> cor_VT;
  std::vector<std::size_t> state;
};
}  // namespace amsim
#endif  // AMSIMCPP_MATING_H
