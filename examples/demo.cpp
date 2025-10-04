#include <iostream>
#include <cstddef>
#include <vector>
#include <string>
#include <amsim/genome.hpp>
#include <amsim/phenotype.hpp>

int main() {
  // Parameters for genome
  std::size_t n_ind = 100000;
  std::size_t n_loc = 5000;
  std::vector<double> v_mut(n_loc, 1e-8);
  std::vector<double> v_rec(n_loc, 0.5);
  std::vector<double> v_maf(n_loc, 0.5);
  int rng_seed = 42;

  std::cout << "Hello from amsim" << std::endl;

  amsim::Genome genome(n_ind, n_loc, v_mut, v_rec, v_maf, rng_seed);

  std::cout << "Initialised genome" << std::endl;

  // Parameters for phenotype
  std::string name = "height";
  std::vector<std::size_t> loci = {1, 2, 3, 4, 5, 6, 7, 8};
  double h2 = 0.5;
  
  amsim::Phenotype phenotype(name, loci, h2);

  std::cout << "Initialised phenotype" << std::endl;

  // Run the simulation
  genome.generate_haplotypes();
  std::cout << "Generated haplotypes" << std::endl;
  genome.transpose();
  std::cout << "Tranposed" << std::endl;
  genome.compute_mafs();
  std::cout << "Computed MAFs" << std::endl;
  genome.compute_stats();
  std::cout << "Computed stats" << std::endl;
  genome.transpose();
  std::cout << "Tranposed" << std::endl;
  phenotype.score(genome);
  std::cout << "Scored phenotype" << std::endl;
}
