#pragma once

#include <amsim/core/params.h>
#include <amsim/data/genome.h>
#include <amsim/data/mating.h>
#include <amsim/data/phenome.h>

namespace amsim {

// State holds the current state of our simulation
struct State {
  std::size_t gen = 0;
  std::size_t parity = 0;

  // store current and parent generation
  std::array<GenoBuf, 2> genos;
  std::array<PhenoBuf, 2> phenos;
  std::array<Matching, 2> matchings;
  std::array<Matching, 2> inv_matchings;

  GenoBuf& geno() { return genos[parity]; }
  PhenoBuf& pheno() { return phenos[parity]; }
  Matching& matching() { return matchings[parity]; }
  Matching& inv_matching() { return inv_matchings[parity]; }

  const GenoBuf& geno() const { return genos[parity]; }
  const PhenoBuf& pheno() const { return phenos[parity]; }
  const Matching& matching() const { return matchings[parity]; }
  const Matching& inv_matching() const { return inv_matchings[parity]; }

  GenoBuf& geno_par() { return genos[1 - parity]; }
  PhenoBuf& pheno_par() { return phenos[1 - parity]; }
  Matching& matching_par() { return matchings[1 - parity]; }
  Matching& inv_matching_par() { return inv_matchings[1 - parity]; }

  const GenoBuf& geno_par() const { return genos[1 - parity]; }
  const PhenoBuf& pheno_par() const { return phenos[1 - parity]; }
  const Matching& matching_par() const { return matchings[1 - parity]; }
  const Matching& inv_matching_par() const {
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

inline State build_state(const Params& params) {
  return State{
    .genos = {GenoBuf(params), GenoBuf(params)},
    .phenos = {PhenoBuf(params), PhenoBuf(params)},
    .matchings = {
        std::vector<std::size_t>(params.geno.n_ind / 2),
        std::vector<std::size_t>(params.geno.n_ind / 2)},
    .inv_matchings = {
        std::vector<std::size_t>(params.geno.n_ind / 2),
        std::vector<std::size_t>(params.geno.n_ind / 2)}
  };
}

}  // namespace amsim
