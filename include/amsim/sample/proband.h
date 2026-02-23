#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>

namespace amsim {

// which generation the probands are in
enum class Generation { Current, Parents };

// proband member sex
enum class Sex : uint8_t { Unknown = 0, Male = 1, Female = 2 };

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
  All = 0b1111,
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
enum class Aggregator { Max, Mean, Min, Identity };

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
       .generation = Generation::Parents,
       .sex = Sex::Male,
       .father = Mate::Unknown,
       .mother = Mate::Unknown},

      {.self = Mate::HusbandMother,
       .generation = Generation::Parents,
       .sex = Sex::Female,
       .father = Mate::Unknown,
       .mother = Mate::Unknown},

      {.self = Mate::WifeFather,
       .generation = Generation::Parents,
       .sex = Sex::Male,
       .father = Mate::Unknown,
       .mother = Mate::Unknown},

      {.self = Mate::WifeMother,
       .generation = Generation::Parents,
       .sex = Sex::Female,
       .father = Mate::Unknown,
       .mother = Mate::Unknown}};
  static constexpr Aggregator AggDefault = Aggregator::Mean;
};

template <>
struct ProbandData<Proband::Family> {
  using ProbandEnum = Family;
  static constexpr std::size_t ProbandSize = 6;
  static constexpr ProbandEnum ProbandDefault = Family::All;
  static constexpr ProbandMember<ProbandEnum> ProbandMembers[ProbandSize] = {
      {.self = Family::Father,
       .generation = Generation::Parents,
       .sex = Sex::Male,
       .father = Family::Unknown,
       .mother = Family::Unknown},

      {.self = Family::Mother,
       .generation = Generation::Parents,
       .sex = Sex::Female,
       .father = Family::Unknown,
       .mother = Family::Unknown},

      {.self = Family::Son,
       .generation = Generation::Current,
       .sex = Sex::Male,
       .father = Family::Father,
       .mother = Family::Mother},

      {.self = Family::SonWife,
       .generation = Generation::Current,
       .sex = Sex::Female,
       .father = Family::Unknown,
       .mother = Family::Unknown},

      {.self = Family::DaughterHusband,
       .generation = Generation::Current,
       .sex = Sex::Male,
       .father = Family::Unknown,
       .mother = Family::Unknown},

      {.self = Family::Daughter,
       .generation = Generation::Current,
       .sex = Sex::Female,
       .father = Family::Father,
       .mother = Family::Mother},
  };
  static constexpr Aggregator AggDefault = Aggregator::Mean;
};

// weighting functions are functions that act on aggregate proband data
// and return
using WeightFunction =
    std::function<void(const Eigen::MatrixXd&, Eigen::VectorXd&)>;

}  // namespace amsim
