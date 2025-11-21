#ifndef AMSIMCPP_PHENOARCH_H
#define AMSIMCPP_PHENOARCH_H

#include <amsim/rng.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace amsim {

class PhenoArch {
 public:
  PhenoArch(
      std::size_t n_pheno,
      std::size_t n_loc_total,
      std::vector<std::size_t> n_loc,
      std::vector<double> h2_gen,
      std::vector<double> h2_env,
      std::vector<double> gen_cor,
      std::vector<double> env_cor,
      const rng::Xoshiro256ss& rng);

  void gen_env(double* ptr_env, std::size_t n_ind);
  void optim_arch(std::size_t max_it, double eps = 1e-12);
  void print_correlations(const std::vector<std::size_t>& intersect) const;
  std::vector<double> env_chol() const noexcept { return env_chol_; }
  std::vector<double> gen_cor() const noexcept { return gen_cor_; }
  std::vector<std::size_t> pheno_mask(std::size_t pheno_id) const;

 private:
  const std::size_t n_pheno_;
  const std::size_t n_loc_tot_;
  const std::vector<std::size_t> n_loc_;
  const std::vector<double> h2_gen_;
  const std::vector<double> h2_env_;
  const std::vector<double> gen_cor_;

  std::vector<double> env_chol_;
  std::vector<std::uint64_t> loc_mask_;

  rng::NormalPolar rng_polar_;
  rng::UniformIntRange rng_unf_;

  std::vector<std::uint64_t> initMask();
  std::vector<double> initEnergy() const;
  std::vector<double> initWeights() const;
  std::vector<std::size_t> initIntersect(
      const std::vector<std::uint64_t>& mask) const;
};

}  // namespace amsim

#endif  // AMSIMCPP_PHENOARCH_H
