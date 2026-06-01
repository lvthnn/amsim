#pragma once

#include <amsim/core/state.h>
#include <amsim/sample/proband.h>

#include <format>

namespace amsim {

template <typename ProbandEnum>
inline std::size_t member_index(
    const State& state, ProbandEnum member, std::size_t proband_id) {
  if constexpr (std::is_same_v<ProbandEnum, Individual>) {
    return proband_id;
  }

  else if constexpr (std::is_same_v<ProbandEnum, Family>) {
    switch (member) {
      case Family::Father:
        return proband_id;
      case Family::Mother:
        return state.matching(Generation::Parents)[proband_id];
      case Family::Son:
        return proband_id;
      case Family::SonWife:
        return state.matching()[proband_id];
      case Family::DaughterHusband:
        return state
            .inv_matching()[state.matching(Generation::Parents)[proband_id]];
      case Family::Daughter:
        return state.n_sex + state.inv_matching()[proband_id];
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
        return state.inv_matching(
            Generation::Parents)[state.matching()[proband_id]];
      case Mate::WifeMother:
        return state.matching(Generation::Parents)[state.inv_matching(
            Generation::Parents)[proband_id]];
      default:
        throw std::invalid_argument("Invalid or compound Mate member type");
    }
  }
}

template <typename ProbandEnum>
inline std::size_t member_index(
    const State& state,
    const ProbandMember<ProbandEnum>& member,
    std::size_t proband_id) {
  return member_index(state, member.self, proband_id);
}

template <typename ProbandEnum>
inline Generation member_generation(ProbandEnum member) {
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
inline Generation member_generation(const ProbandMember<ProbandEnum>& member) {
  return member.generation;
}

template <typename ProbandEnum>
inline std::uint8_t member_geno_plink(
    const State& state,
    ProbandEnum member,
    std::size_t proband_id,
    std::size_t locus) {
  Generation generation = member_generation(member);
  std::size_t index = member_index(state, member, proband_id);

  std::size_t word = index / 64;
  std::size_t bit = index % 64;

  std::uint8_t h0 = (state.geno(generation).h0()(locus, word) >> bit) & 1;
  std::uint8_t h1 = (state.geno(generation).h1()(locus, word) >> bit) & 1;
  return (h0 & h1) | ((h0 | h1) << 1);
}

template <typename ProbandEnum>
inline std::uint8_t member_geno_plink(
    const State& state,
    const ProbandMember<ProbandEnum>& member,
    std::size_t proband_id,
    std::size_t locus) {
  std::size_t index = member_index(state, member, proband_id);

  std::size_t word = index / 64;
  std::size_t bit = index % 64;

  std::uint8_t h0 =
      (state.geno(member.generation).h0()(locus, word) >> bit) & 1;
  std::uint8_t h1 =
      (state.geno(member.generation).h1()(locus, word) >> bit) & 1;
  return (h0 & h1) | ((h0 | h1) << 1);
}

template <typename ProbandEnum>
inline double member_pheno(
    const State& state,
    ProbandEnum member,
    std::size_t proband_id,
    std::size_t pheno_id,
    Component component = Component::Total) {
  Generation generation = member_generation(member);
  std::size_t index = member_index(state, member, proband_id);
  return state.pheno(generation)(pheno_id, component)(index);
}

template <typename ProbandEnum>
inline double member_pheno(
    const State& state,
    const ProbandMember<ProbandEnum>& member,
    std::size_t proband_id,
    std::size_t pheno_id,
    Component component = Component::Total) {
  std::size_t index = member_index(state, member, proband_id);
  return state.pheno(member.generation)(pheno_id, component)(index);
}

template <Proband P>
constexpr std::string_view prefix_of(
    typename ProbandData<P>::ProbandEnum member);

template <>
constexpr std::string_view prefix_of<Proband::Individual>(Individual member) {
  switch (member) {
    case Individual::Self:
      return "IND";
    default:
      return "";
  }
}

template <>
constexpr std::string_view prefix_of<Proband::Mate>(Mate member) {
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
constexpr std::string_view prefix_of<Proband::Family>(Family member) {
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
inline std::string member_id(
    const State& state,
    typename ProbandData<P>::ProbandEnum member,
    std::size_t proband_id) {
  if (member == ProbandData<P>::ProbandEnum::Unknown) return "0";
  return std::format(
      "{}{}", prefix_of<P>(member), member_index(state, member, proband_id));
}

}  // namespace amsim
