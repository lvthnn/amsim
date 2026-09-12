// This file is part of amsim, copyright (C) 2025-2026 Kári Hlynsson.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <amsim/core/generation.h>
#include <amsim/core/params.h>
#include <amsim/data/genome.h>
#include <amsim/data/mating.h>
#include <amsim/data/pedigree.h>
#include <amsim/data/phenome.h>

namespace amsim {

// State holds the current state of our simulation
struct State {
  const std::size_t n_sex;
  std::size_t gen = 0;
  std::size_t rep = 0;
  std::size_t parity = 0;

  // store current and parent generation
  std::array<GenoBuf, 2> genos;
  std::array<PhenoBuf, 2> phenos;
  std::array<Matching, 2> matchings;
  std::array<Matching, 2> inv_matchings;

  // for stoing family data
  Pedigree pedigree;

  std::size_t getParity(Generation generation) const {
    return (generation == Generation::Current) ? parity : 1 - parity;
  }

  GenoBuf& geno(Generation generation = Generation::Current) {
    return genos[getParity(generation)];
  }

  PhenoBuf& pheno(Generation generation = Generation::Current) {
    return phenos[getParity(generation)];
  }

  Matching& matching(Generation generation = Generation::Current) {
    return matchings[getParity(generation)];
  }

  Matching& invMatching(Generation generation = Generation::Current) {
    return inv_matchings[getParity(generation)];
  }

  const GenoBuf& geno(Generation generation = Generation::Current) const {
    return genos[getParity(generation)];
  }

  const PhenoBuf& pheno(Generation generation = Generation::Current) const {
    return phenos[getParity(generation)];
  }

  const Matching& matching(Generation generation = Generation::Current) const {
    return matchings[getParity(generation)];
  }

  const Matching& invMatching(
      Generation generation = Generation::Current) const {
    return inv_matchings[getParity(generation)];
  }

  void transpose() {
    (*this).geno(Generation::Current).transpose();
    (*this).geno(Generation::Parents).transpose();
  }

  void advance() {
    gen += 1;
    parity = gen % 2;
  }

  void updatePedigree() { pedigree.push(matching(), invMatching()); }
};

inline State buildState(const Params& params) {
  return State{
      .n_sex = params.global.n_ind / 2,
      .rep = params.global.rep_id,
      .genos = {GenoBuf(params), GenoBuf(params)},
      .phenos = {PhenoBuf(params), PhenoBuf(params)},
      .matchings =
          {std::vector<std::size_t>(params.global.n_ind / 2),
           std::vector<std::size_t>(params.global.n_ind / 2)},
      .inv_matchings =
          {std::vector<std::size_t>(params.global.n_ind / 2),
           std::vector<std::size_t>(params.global.n_ind / 2)},
      .pedigree = Pedigree(params)};
}

}  // namespace amsim
