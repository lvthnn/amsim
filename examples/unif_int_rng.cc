//
// Created by Kári Hlynsson on 8.10.2025.
//

#include <cstddef>
#include <iostream>
#include <vector>

#include <amsim/rng.h>

int main() {
  constexpr std::size_t n_samples = 1e6;
  constexpr std::size_t range_upper = 50;
  amsim::rng::Xoshiro256ss rng = amsim::rng::seed_xoshiro(amsim::rng::auto_seed(12345ull));
  amsim::rng::UniformIntRange unif_dist(rng);

  for (std::size_t i = 0; i < n_samples; i++)
    std::cout << unif_dist.sample(range_upper) << "\n";
}