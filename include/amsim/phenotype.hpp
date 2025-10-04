#pragma once
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include <amsim/genome.hpp>

namespace amsim {
  struct Phenotype {
    const std::string name;
    const std::vector<size_t> loci;
    const std::vector<uint64_t> word_masks;
    const double h2;

    std::vector<double> loc_effects;
    std::unordered_map<std::size_t, std::uint64_t> loc_mask;

    std::vector<double> values_gen;
    std::vector<double> values_env;
    std::vector<double> values;

		Phenotype(const std::string name_, const std::vector<size_t>& loci_,
              const double h2_);

    void score(Genome& genome);
  };
}
