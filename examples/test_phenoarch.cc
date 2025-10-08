//
// Created by Kári Hlynsson on 7.10.2025.
//

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
  std::size_t seed = 123456ull;
  std::vector<double> h2_gen{0.5, 0.5, 0.5};
  std::vector<std::size_t> n_loc{3000, 3000, 3000};
  std::vector<double> env_buffer(n_pheno * n_ind);

  // specify column-major layout
  std::vector<double> gen_cor{
    1.0, 0.2, 0.1,   // col 1
    0.2, 1.0, 0.1,   // col 2
    0.1, 0.1, 1.0    // col 3
  };

  std::vector<double> env_cor{
    1.0, 0.2, 0.1,
    0.2, 1.0, 0.1,
    0.1, 0.1, 1.0
  };

  amsim::rng::Xoshiro256ss rng = amsim::rng::seed_xoshiro(amsim::rng::auto_seed(seed));

  amsim::PhenoArch arch(n_pheno, n_loc, 9000, h2_gen, gen_cor,
                      env_cor, rng);

  arch.optim_arch();

  // std::cout << std::setprecision(5);
  // std::size_t cnt = 0;
  // for (const double val : env_buffer) {
  //   std::cout << val << "\t";
  //   cnt = (cnt + 1) % n_pheno;
  //   if (cnt == 0) std::cout << "\n";
  // }
}