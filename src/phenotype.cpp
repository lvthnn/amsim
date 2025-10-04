#include <amsim/phenotype.hpp>
#include <amsim/genome.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <unordered_map>

namespace amsim {
  
  Phenotype::Phenotype(const std::string name_,
                       const std::vector<std::size_t>& loci_,
                       const double h2_)
    : name(name_),
      loci(std::move(loci_)),
      h2(h2_) {
    for (std::size_t loc : loci_) {
      std::size_t block = loc / 64;
      std::size_t offset = loc % 64;
      loc_mask[block] |= std::uint64_t(1) << offset;
    }

    loc_effects.resize(loci.size(), std::sqrt(h2 / loci.size()));
  }

  void Phenotype::score(Genome& genome) {
    std::size_t n_ind = genome.H0().n_ind();
    HapMat& H0 = genome.H0();
    HapMat& H1 = genome.H1();

    values_gen.assign(n_ind, 0.0);

    for (std::size_t el = 0; el < loci.size(); el++) {
      std::size_t loc = loci[el];
      std::size_t block = loc / 64;
      std::size_t offset = loc % 64;\
      double loc_sd = std::sqrt(genome.v_lvar(loc));
      double loc_mean = genome.v_lmean(loc);

      for (std::size_t ind = 0; ind < n_ind; ind++) {
        int h0_loc = (H0(ind, block) >> offset) & 1ull;
        int h1_loc = (H1(ind, block) >> offset) & 1ull;
        int geno_loc = h0_loc + h1_loc;

        values_gen[ind] += loc_effects[loc] * (geno_loc - loc_mean) / loc_sd;
      }
    }
  }
}
