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

#include <amsim/core/state.h>
#include <amsim/sample/proband.h>

#include <format>

namespace amsim {

template <typename ProbandEnum>
inline std::size_t memberIndex(
    const State& state, ProbandEnum member, std::size_t proband_id) {
  if constexpr (std::is_same_v<ProbandEnum, Individual>) {
    return proband_id;
  }

  else if constexpr (std::is_same_v<ProbandEnum, Family>) {
    switch (member) {
      case Family::Father:
        return proband_id;
      case Family::Mother:
        return state.n_sex + state.matching(Generation::Parents)[proband_id];
      case Family::Son:
        return proband_id;
      case Family::SonWife:
        return state.n_sex + state.matching()[proband_id];
      case Family::DaughterHusband:
        return state
            .invMatching()[state.matching(Generation::Parents)[proband_id]];
      case Family::Daughter:
        return state.n_sex + state.matching(Generation::Parents)[proband_id];
      default:
        throw std::invalid_argument("Invalid or compound Family member type");
    }
  }

  else if constexpr (std::is_same_v<ProbandEnum, Mate>) {
    switch (member) {
      case Mate::Husband:
        return proband_id;
      case Mate::Wife:
        return state.n_sex + state.matching()[proband_id];
      case Mate::HusbandFather:
        return proband_id;
      case Mate::HusbandMother:
        return state.matching(Generation::Parents)[proband_id];
      case Mate::WifeFather:
        return state.invMatching(
            Generation::Parents)[state.matching()[proband_id]];
      case Mate::WifeMother:
        return state.matching(Generation::Parents)[state.invMatching(
            Generation::Parents)[proband_id]];
      default:
        throw std::invalid_argument("Invalid or compound Mate member type");
    }
  }
}

template <typename ProbandEnum>
inline std::size_t memberIndex(
    const State& state,
    const ProbandMember<ProbandEnum>& member,
    std::size_t proband_id) {
  return memberIndex(state, member.self, proband_id);
}

template <typename ProbandEnum>
inline Generation memberGeneration(ProbandEnum member) {
  if constexpr (std::is_same_v<ProbandEnum, Individual>) {
    return Generation::Current;
  } else if constexpr (std::is_same_v<ProbandEnum, Family>) {
    switch (member) {
      case Family::Father:
      case Family::Mother:
        return Generation::Parents;
      default:
        return Generation::Current;
    }
  } else if constexpr (std::is_same_v<ProbandEnum, Mate>) {
    switch (member) {
      case Mate::Husband:
      case Mate::Wife:
        return Generation::Current;
      default:
        return Generation::Parents;
    }
  }
}

template <typename ProbandEnum>
inline Generation memberGeneration(const ProbandMember<ProbandEnum>& member) {
  return member.generation;
}

template <typename ProbandEnum>
inline std::uint8_t memberGenoPLINK(
    const State& state,
    ProbandEnum member,
    std::size_t proband_id,
    std::size_t locus) {
  Generation generation = memberGeneration(member);
  std::size_t index = memberIndex(state, member, proband_id);

  std::size_t word = index / 64;
  std::size_t bit = index % 64;

  std::uint8_t h0 = (state.geno(generation).h0()(locus, word) >> bit) & 1;
  std::uint8_t h1 = (state.geno(generation).h1()(locus, word) >> bit) & 1;
  return (h0 & h1) | ((h0 | h1) << 1);
}

template <typename ProbandEnum>
inline std::uint8_t memberGenoPLINK(
    const State& state,
    const ProbandMember<ProbandEnum>& member,
    std::size_t proband_id,
    std::size_t locus) {
  std::size_t index = memberIndex(state, member, proband_id);

  std::size_t word = index / 64;
  std::size_t bit = index % 64;

  std::uint8_t h0 =
      (state.geno(member.generation).h0()(locus, word) >> bit) & 1;
  std::uint8_t h1 =
      (state.geno(member.generation).h1()(locus, word) >> bit) & 1;
  return (h0 & h1) | ((h0 | h1) << 1);
}

template <typename ProbandEnum>
inline double memberPheno(
    const State& state,
    ProbandEnum member,
    std::size_t proband_id,
    std::size_t pheno_id,
    Component component = Component::Total) {
  Generation generation = memberGeneration(member);
  std::size_t index = memberIndex(state, member, proband_id);
  return state.pheno(generation)(pheno_id, component)(index);
}

template <typename ProbandEnum>
inline double memberPheno(
    const State& state,
    const ProbandMember<ProbandEnum>& member,
    std::size_t proband_id,
    std::size_t pheno_id,
    Component component = Component::Total) {
  std::size_t index = memberIndex(state, member, proband_id);
  return state.pheno(member.generation)(pheno_id, component)(index);
}

template <Proband P>
constexpr std::string_view prefixOf(
    typename ProbandData<P>::ProbandEnum member);

template <>
constexpr std::string_view prefixOf<Proband::Individual>(Individual member) {
  switch (member) {
    case Individual::Self:
      return "IND";
    default:
      return "";
  }
}

template <>
constexpr std::string_view prefixOf<Proband::Mate>(Mate member) {
  switch (member) {
    case Mate::Husband:
      return "HUS";
    case Mate::Wife:
      return "WIF";
    case Mate::HusbandFather:
      return "HUF";
    case Mate::HusbandMother:
      return "HUM";
    case Mate::WifeFather:
      return "WFF";
    case Mate::WifeMother:
      return "WFM";
    default:
      return "";
  }
}

template <>
constexpr std::string_view prefixOf<Proband::Family>(Family member) {
  switch (member) {
    case Family::Father:
      return "FAT";
    case Family::Mother:
      return "MOT";
    case Family::Son:
      return "SON";
    case Family::Daughter:
      return "DAU";
    case Family::SonWife:
      return "SOW";
    case Family::DaughterHusband:
      return "DAH";
    default:
      return "";
  }
}

template <Proband P>
inline std::string memberID(
    const State& state,
    typename ProbandData<P>::ProbandEnum member,
    std::size_t proband_id) {
  if (member == ProbandData<P>::ProbandEnum::Unknown) return "0";
  return std::format(
      "{}{}", prefixOf<P>(member), memberIndex(state, member, proband_id));
}

}  // namespace amsim
