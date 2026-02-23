#pragma once

#include <amsim/core/params.h>
#include <amsim/data/genome.h>
#include <amsim/data/mating.h>
#include <amsim/data/phenome.h>
#include <amsim/sample/proband.h>

namespace amsim {

// State holds the current state of our simulation
struct State {
  const std::size_t n_sex;
  std::size_t gen = 0;
  std::size_t parity = 0;

  // store current and parent generation
  std::array<GenoBuf, 2> genos;
  std::array<PhenoBuf, 2> phenos;
  std::array<Matching, 2> matchings;
  std::array<Matching, 2> inv_matchings;

  std::size_t get_parity(Generation generation) const {
    return (generation == Generation::Current) ? parity : 1 - parity;
  }

  GenoBuf& geno(Generation generation = Generation::Current) {
    return genos[get_parity(generation)];
  }

  PhenoBuf& pheno(Generation generation = Generation::Current) {
    return phenos[get_parity(generation)];
  }

  Matching& matching(Generation generation = Generation::Current) {
    return matchings[get_parity(generation)];
  }

  Matching& inv_matching(Generation generation = Generation::Current) {
    return inv_matchings[get_parity(generation)];
  }

  const GenoBuf& geno(Generation generation = Generation::Current) const {
    return genos[get_parity(generation)];
  }

  const PhenoBuf& pheno(Generation generation = Generation::Current) const {
    return phenos[get_parity(generation)];
  }

  const Matching& matching(Generation generation = Generation::Current) const {
    return matchings[get_parity(generation)];
  }

  const Matching& inv_matching(
      Generation generation = Generation::Current) const {
    return inv_matchings[get_parity(generation)];
  }

  void transpose() {
    (*this).geno(Generation::Current).transpose();
    (*this).geno(Generation::Parents).transpose();
  }

  void advance() {
    gen += 1;
    parity = gen % 2;
  }
};

inline State build_state(const Params& params) {
  return State{
      .n_sex = params.geno.n_ind / 2,
      .genos = {GenoBuf(params), GenoBuf(params)},
      .phenos = {PhenoBuf(params), PhenoBuf(params)},
      .matchings =
          {std::vector<std::size_t>(params.geno.n_ind / 2),
           std::vector<std::size_t>(params.geno.n_ind / 2)},
      .inv_matchings = {
          std::vector<std::size_t>(params.geno.n_ind / 2),
          std::vector<std::size_t>(params.geno.n_ind / 2)}};
}

}  // namespace amsim
