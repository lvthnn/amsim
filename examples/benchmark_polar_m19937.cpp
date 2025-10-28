//
// Created by Kári Hlynsson on 7.10.2025.
//

#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>

#include <amsim/rng.h>
#include <amsim/utils.h>

int main(int argc, char *argv[]) {
  std::size_t n_draws = 256000 * 50;
  std::random_device rd;
  std::mt19937 g(rd());
  std::normal_distribution<double> n(0, 1);
  std::vector<double> rng_vec(n_draws);

  amsim::utils::time_step("Standard library RNG", [&]{
    for (std::size_t i = 0; i < n_draws; i++)
      rng_vec[i] = n(g);
  });

  amsim::rng::Xoshiro256ss rng = amsim::rng::seed_xoshiro(amsim::rng::auto_seed(12345ull));
  amsim::rng::NormalPolar polar(rng);

  amsim::utils::time_step("amsim Polar RNG", [&]{ polar.fill(rng_vec.data(), rng_vec.size()); });
}