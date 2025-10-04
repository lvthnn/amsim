#include <iostream>
#include <cstddef>
#include <vector>
#include <string>
#include <numeric>
#include <chrono>

#include <amsim/genome.hpp>
#include <amsim/phenotype.hpp>

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

double variance(std::vector<double> a, int n)
{
  // Compute mean (average of elements)
  double sum = 0;
  for (int i = 0; i < n; i++)
      sum += a[i];
  double mean = (double)sum / (double)n;

  // Compute sum squared 
  // differences with mean.
  double sqDiff = 0;
  for (int i = 0; i < n; i++) 
      sqDiff += (a[i] - mean) * 
                (a[i] - mean);
  return sqDiff / n;
}

std::vector<std::size_t> range_inclusive_iota(std::size_t a, std::size_t b) {
  std::vector<std::size_t> v(b - a + 1);
  std::iota(v.begin(), v.end(), a);  // fills a, a+1, ..., b
  return v;
}

int main() {
  // Parameters for genome
  std::size_t n_ind = 256000;
  std::size_t n_loc = 14000;
  std::vector<double> v_mut(n_loc, 1e-8);
  std::vector<double> v_rec(n_loc, 0.5);
  std::vector<double> v_maf(n_loc, 0.5);
  int rng_seed = 42;

  std::cout << "Hello from amsim" << std::endl;

  amsim::Genome genome(n_ind, n_loc, v_mut, v_rec, v_maf, rng_seed);

  std::cout << "Initialised genome" << std::endl;

  // Parameters for phenotype
  std::string name = "height";
  std::vector<std::size_t> loci = range_inclusive_iota(0, 999);
  double h2 = 0.5;

  amsim::Phenotype phenotype(name, loci, h2);

  std::cout << "Initialised phenotype" << std::endl;

  // Run the simulation
  // Redo this section, adding timers around each step
  const auto total_start = SteadyClock::now();

  time_step("Generated haplotypes", [&]{ genome.generate_haplotypes(); });
  time_step("Transposed",          [&]{ genome.transpose(); });
  time_step("Computed MAFs",       [&]{ genome.compute_mafs(); });
  time_step("Computed stats",      [&]{ genome.compute_stats(); });
  time_step("Transposed",          [&]{ genome.transpose(); });
  time_step("Scored phenotype",    [&]{ phenotype.score(genome); });

  const auto total_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          SteadyClock::now() - total_start).count();
  std::cout << "Total elapsed: " << total_ms << " ms\n\n";
}
