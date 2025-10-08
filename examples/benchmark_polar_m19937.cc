//
// Created by Kári Hlynsson on 7.10.2025.
//

#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>

#include <amsim/rng.h>

using SteadyClock = std::chrono::steady_clock;

template <class F>
void time_step(const char* label, F&& fn) {
  const auto t0 = SteadyClock::now();
  fn();
  const auto t1 = SteadyClock::now();
  const auto ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  std::cout << label << " (" << ms << " ms)\n";
}

int main(int argc, char *argv[]) {
  std::size_t n_draws = 256000 * 50;
  std::random_device rd;
  std::mt19937 g(rd());
  std::normal_distribution<double> n(0, 1);
  std::vector<double> rng_vec(n_draws);

  time_step("Standard library RNG", [&]{
    for (std::size_t i = 0; i < n_draws; i++)
      rng_vec[i] = n(g);
  });

  amsim::rng::Xoshiro256ss rng = amsim::rng::seed_xoshiro(amsim::rng::auto_seed(12345ull));
  amsim::rng::NormalPolar polar(rng);

  time_step("amsim Polar RNG", [&]{ polar.fill(rng_vec.data(), rng_vec.size()); });
}