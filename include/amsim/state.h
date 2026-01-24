#pragma once

#include <amsim/data/genome.h>
#include <amsim/data/phenome.h>
#include <amsim/data/mating.h>

namespace amsim {

// State holds the current state of our simulation
struct State {
  genome::GenoBuf& geno;
  phenome::PhenoBuf& pheno;
  mating::Matching& matching;
};

}
