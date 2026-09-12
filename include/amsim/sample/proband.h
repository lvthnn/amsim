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
#include <amsim/core/utils.h>
#include <amsim/sample/weight.h>

#include <boost/algorithm/string/case_conv.hpp>

namespace amsim {

// proband types — these define the unit of sampling
enum class Proband { Individual, Mate, Family };

// proband subtypes — the members constituting a sampling unit and various
// combinations of them, such as parents, siblings, in-laws, etc.
enum class Individual : uint8_t { Unknown = 0b0, Self = 0b1 };

enum class Mate : uint8_t {
  Unknown = 0b000000,
  Husband = 0b000001,
  Wife = 0b000010,
  HusbandFather = 0b000100,
  HusbandMother = 0b001000,
  WifeFather = 0b010000,
  WifeMother = 0b100000,
  All = 0b111111,
  Couple = Husband | Wife,
  HusbandInLaws = WifeFather | WifeMother,
  WifeInLaws = HusbandFather | HusbandMother,
  Parents = HusbandInLaws | WifeInLaws,
  HusbandFamily = Husband | WifeInLaws,
  WifeFamily = Wife | HusbandInLaws,
  Males = Husband | HusbandFather | WifeFather,
  Females = Wife | HusbandMother | WifeMother
};

enum class Family : uint8_t {
  Unknown = 0b000000,
  Father = 0b000001,
  Mother = 0b000010,
  Son = 0b000100,
  Daughter = 0b001000,
  SonWife = 0b010000,
  DaughterHusband = 0b100000,
  All = 0b111111,
  Parents = Father | Mother,
  Siblings = Son | Daughter,
  Males = Father | Son,
  Females = Mother | Daughter
};

constexpr uint8_t operator&(Individual a, Individual b) {
  return static_cast<uint8_t>(a) & static_cast<uint8_t>(b);
}

constexpr uint8_t operator&(Mate a, Mate b) {
  return static_cast<uint8_t>(a) & static_cast<uint8_t>(b);
}

constexpr uint8_t operator&(Family a, Family b) {
  return static_cast<uint8_t>(a) & static_cast<uint8_t>(b);
}

constexpr Individual operator|(Individual a, Individual b) {
  return static_cast<Individual>(
      static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
}

constexpr Mate operator|(Mate a, Mate b) {
  return static_cast<Mate>(static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
}

constexpr Family operator|(Family a, Family b) {
  return static_cast<Family>(static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
}

// aggregator types — functions to reduce a proband subtype into a statistics
// that can enter into a sampling probability transformer
enum class Aggregator { Max, Min, Mean, Identity };

inline Aggregator aggregatorFromString(const std::string& s) {
  std::string l = boost::to_lower_copy(s);
  if (l == "max") return Aggregator::Max;
  if (l == "min") return Aggregator::Min;
  if (l == "mean") return Aggregator::Mean;
  if (l == "identity") return Aggregator::Identity;
  throw std::runtime_error("Unknown aggregator type" + l);
}

// declares member information
template <typename ProbandEnum>
struct ProbandMember {
  ProbandEnum self;
  ProbandEnum father;
  ProbandEnum mother;
  Generation generation;
  Sex sex;
};

// for each kind of proband type, the extracted data from sampling has a
// fixed layout, depending on the size (or "dimension") of the sampling unit
// with respect to the number of "atomic" members comprising it
template <Proband P>
struct ProbandData;

template <>
struct ProbandData<Proband::Individual> {
  using ProbandEnum = Individual;
  static constexpr std::size_t ProbandSize = 1;
  static constexpr ProbandEnum ProbandDefault = Individual::Self;
  static constexpr ProbandMember<ProbandEnum> ProbandMembers[ProbandSize] = {
      {.self = Individual::Self,
       .generation = Generation::Current,
       .sex = Sex::Unknown}};
  static constexpr Aggregator AggDefault = Aggregator::Identity;
  static ProbandEnum fromString(const std::string& s) {
    std::string l = boost::to_lower_copy(s);
    if (l == "self") return Individual::Self;
    throw std::runtime_error("Unknown Individual proband member " + s);
  }
};

template <>
struct ProbandData<Proband::Mate> {
  using ProbandEnum = Mate;
  static constexpr std::size_t ProbandSize = 6;
  static constexpr ProbandEnum ProbandDefault = Mate::Couple;
  static constexpr ProbandMember<ProbandEnum> ProbandMembers[ProbandSize] = {
      {.self = Mate::Husband,
       .father = Mate::HusbandFather,
       .mother = Mate::HusbandMother,
       .generation = Generation::Current,
       .sex = Sex::Male},

      {.self = Mate::Wife,
       .father = Mate::WifeFather,
       .mother = Mate::WifeMother,
       .generation = Generation::Current,
       .sex = Sex::Female},

      {.self = Mate::HusbandFather,
       .father = Mate::Unknown,
       .mother = Mate::Unknown,
       .generation = Generation::Parents,
       .sex = Sex::Male},

      {.self = Mate::HusbandMother,
       .father = Mate::Unknown,
       .mother = Mate::Unknown,
       .generation = Generation::Parents,
       .sex = Sex::Female},

      {.self = Mate::WifeFather,
       .father = Mate::Unknown,
       .mother = Mate::Unknown,
       .generation = Generation::Parents,
       .sex = Sex::Male},

      {.self = Mate::WifeMother,
       .father = Mate::Unknown,
       .mother = Mate::Unknown,
       .generation = Generation::Parents,
       .sex = Sex::Female}};
  static constexpr Aggregator AggDefault = Aggregator::Mean;
  static ProbandEnum fromString(const std::string& s) {
    std::string l = boost::to_lower_copy(s);
    if (l == "husband") return Mate::Husband;
    if (l == "wife") return Mate::Wife;
    if (l == "husbandfather") return Mate::HusbandFather;
    if (l == "husbandmother") return Mate::HusbandMother;
    if (l == "wifefather") return Mate::WifeFather;
    if (l == "wifemother") return Mate::WifeMother;
    if (l == "all") return Mate::All;
    if (l == "couple") return Mate::Couple;
    if (l == "husbandinlaws") return Mate::HusbandInLaws;
    if (l == "wifeinlaws") return Mate::WifeInLaws;
    if (l == "parents") return Mate::Parents;
    if (l == "husbandfamily") return Mate::HusbandFamily;
    if (l == "wifefamily") return Mate::WifeFamily;
    if (l == "males") return Mate::Males;
    if (l == "females") return Mate::Females;
    throw std::runtime_error("Unknown Mate proband member " + l);
  }
};

template <>
struct ProbandData<Proband::Family> {
  using ProbandEnum = Family;
  static constexpr std::size_t ProbandSize = 6;
  static constexpr ProbandEnum ProbandDefault = Family::All;
  static constexpr ProbandMember<ProbandEnum> ProbandMembers[ProbandSize] = {
      {.self = Family::Father,
       .father = Family::Unknown,
       .mother = Family::Unknown,
       .generation = Generation::Parents,
       .sex = Sex::Male},

      {.self = Family::Mother,
       .father = Family::Unknown,
       .mother = Family::Unknown,
       .generation = Generation::Parents,
       .sex = Sex::Female},

      {.self = Family::Son,
       .father = Family::Father,
       .mother = Family::Mother,
       .generation = Generation::Current,
       .sex = Sex::Male},

      {.self = Family::SonWife,
       .father = Family::Unknown,
       .mother = Family::Unknown,
       .generation = Generation::Current,
       .sex = Sex::Female},

      {.self = Family::DaughterHusband,
       .father = Family::Unknown,
       .mother = Family::Unknown,
       .generation = Generation::Current,
       .sex = Sex::Male},

      {.self = Family::Daughter,
       .father = Family::Father,
       .mother = Family::Mother,
       .generation = Generation::Current,
       .sex = Sex::Female},
  };
  static constexpr Aggregator AggDefault = Aggregator::Mean;
  static ProbandEnum fromString(const std::string& s) {
    std::string l = boost::to_lower_copy(s);
    if (l == "father") return Family::Father;
    if (l == "mother") return Family::Mother;
    if (l == "son") return Family::Son;
    if (l == "daughter") return Family::Daughter;
    if (l == "sonwife") return Family::SonWife;
    if (l == "daughterhusband") return Family::DaughterHusband;
    if (l == "all") return Family::All;
    if (l == "parents") return Family::Parents;
    if (l == "siblings") return Family::Siblings;
    if (l == "males") return Family::Males;
    if (l == "females") return Family::Females;
    throw std::runtime_error("Unknown Family proband member " + l);
  }
};

template <Proband P>
typename ProbandData<P>::ProbandEnum parseProband(
    const std::vector<std::string>& probands_str) {
  using ProbandEnum = typename ProbandData<P>::ProbandEnum;
  std::vector<ProbandEnum> probands(probands_str.size());

  if (probands.empty())
    throw std::runtime_error("Empty proband string supplied");

  std::ranges::transform(
      probands_str, probands.begin(), [](const std::string& s) {
        return ProbandData<P>::fromString(s);
      });

  return std::accumulate(
      probands.begin() + 1,
      probands.end(),
      probands[0],
      [](const ProbandEnum& a, ProbandEnum b) { return a | b; });
}

}  // namespace amsim
