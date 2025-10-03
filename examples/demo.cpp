#include <iostream>
#include <amsim/genome.hpp>

int main() {
  amsim::Genome g(1000, 5000, {}, {}, {}, 123);
  g.generate_haplotypes();
  std::cout << "yep" << std::endl;
}
