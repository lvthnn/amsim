//------------------------------------------------------------------------------
// amsimcpp : phenoarch.h
//------------------------------------------------------------------------------

#ifndef AMSIMCPP_PHENOARCH_H
#define AMSIMCPP_PHENOARCH_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include <amsim/rng.h>

namespace amsim {

  class PhenoArch {
  public:
    PhenoArch(std::size_t n_pheno, std::vector<std::size_t> n_loc,
              std::size_t n_loc_total, std::vector<double> h2_gen,
              std::vector<double> gen_cor, std::vector<double> env_cor,
              const rng::Xoshiro256ss &rng);
    void gen_env(double* ptr_env, std::size_t n_ind);
    void optim_arch(double eps = 1e-12, std::size_t max_it = 1e5);

    std::vector<double> env_chol() const noexcept { return env_chol_; }
    std::vector<double> gen_cor() const noexcept { return gen_cor_; }

  private:
    const std::size_t n_pheno_;
    const std::vector<std::size_t> n_loc_;
    const std::size_t n_loc_tot_;
    const std::vector<double> h2_gen_;
    const std::vector<double> gen_cor_;

    std::vector<double> env_chol_;
    std::vector<uint64_t> loc_mask_;

    rng::NormalPolar rng_polar_;
    rng::UniformIntRange rng_unf_;
  };

}

#endif //AMSIMCPP_PHENOARCH_H