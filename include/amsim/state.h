#pragma once

#include <amsim/data/genome.h>
#include <amsim/data/mating.h>
#include <amsim/data/phenome.h>
#include <amsim/params.h>

namespace amsim {

// State holds the current state of our simulation
struct State {
  std::size_t gen = 0;
  std::size_t parity = 0;

  // store current and parent generation
  std::array<genome::GenoBuf, 2> genos;
  std::array<phenome::PhenoBuf, 2> phenos;
  std::array<mating::Matching, 2> matchings;
  std::array<mating::Matching, 2> inv_matchings;

  genome::GenoBuf& geno() { return genos[parity]; }
  phenome::PhenoBuf& pheno() { return phenos[parity]; }
  mating::Matching& matching() { return matchings[parity]; }
  mating::Matching& inv_matching() { return inv_matchings[parity]; }

  const genome::GenoBuf& geno() const { return genos[parity]; }
  const phenome::PhenoBuf& pheno() const { return phenos[parity]; }
  const mating::Matching& matching() const { return matchings[parity]; }
  const mating::Matching& inv_matching() const { return inv_matchings[parity]; }

  genome::GenoBuf& geno_par() { return genos[1 - parity]; }
  phenome::PhenoBuf& pheno_par() { return phenos[1 - parity]; }
  mating::Matching& matching_par() { return matchings[1 - parity]; }
  mating::Matching& inv_matching_par() { return inv_matchings[1 - parity]; }

  const genome::GenoBuf& geno_par() const { return genos[1 - parity]; }
  const phenome::PhenoBuf& pheno_par() const { return phenos[1 - parity]; }
  const mating::Matching& matching_par() const { return matchings[1 - parity]; }
  const mating::Matching& inv_matching_par() const {
    return matchings[1 - parity];
  }

  void transpose() {
    (*this).geno().transpose();
    (*this).geno_par().transpose();
  }

  void advance() {
    gen += 1;
    parity = gen % 2;
  }
};

State build_state(const Params& params);

}  // namespace amsim
