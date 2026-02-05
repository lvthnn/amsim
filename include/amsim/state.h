#pragma once

#include <amsim/data/genome.h>
#include <amsim/data/mating.h>
#include <amsim/data/phenome.h>

namespace amsim {

// State holds the current state of our simulation
struct State {
  std::size_t gen = 0;
  genome::GenoBuf& geno;
  phenome::PhenoBuf& pheno;
  mating::Matching& matching;
};

}  // namespace amsim
