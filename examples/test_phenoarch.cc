#include <cstddef>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>

#include <Accelerate/Accelerate.h>

#include <amsim/phenoarch.h>

int main() {
  std::size_t n_pheno = 3;
  std::size_t n_ind = 256000;
  std::size_t seed = 12345678ull;
  std::vector<double> h2_gen{0.5, 0.5, 0.5};
  std::vector<double> h2_env{0.5, 0.5, 0.5};
  std::vector<std::size_t> n_loc{1000, 1000, 1000};
  std::vector<double> env_buffer(n_pheno * n_ind);

  // specify column-major layout
  std::vector<double> gen_cor{
    1.0, 0.2, 0.3,   // col 1
    0.2, 1.0, 0.4,   // col 2
    0.3, 0.4, 1.0    // col 3
  };

  std::vector<double> env_cor{
    1.0, 0.2, 0.1,   // col 1
    0.2, 1.0, 0.1,   // col 2
    0.1, 0.1, 1.0    // col 3
  };

  const amsim::rng::Xoshiro256ss rng = amsim::rng::seed_xoshiro(amsim::rng::auto_seed(seed));

  amsim::PhenoArch arch(n_pheno, 3000, n_loc, h2_gen, h2_env, gen_cor, env_cor, rng);

  arch.optim_arch(std::size_t(3e3));
}
