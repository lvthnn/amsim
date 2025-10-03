#pragma once
#include <cstddef>
#include <string>
#include <vector>
#include <cmath>

namespace amsim {
  struct Phenotype {
    const std::string name;
    const std::vector<size_t> loci;
    const std::vector<uint64_t> word_masks;
    const double h2;
    double scale_const; 

    std::vector<double> values_gen;
    std::vector<double> values_env;

		Phenotype(const std::string name_, const std::vector<size_t>& loci_,
              const double h2_)
    : name(name_),
      loci(loci_),
      h2(h2_) { };

    inline std::vector<double> values() {
      std::vector<double> values(values_gen.size());
      for (std::size_t el = 0; el < values_gen.size(); el++)
        values[el] = values_gen[el] + values_env[el];
      return values;
    }
  };
}
