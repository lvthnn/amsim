#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>

#include <vector>

namespace amsim {

class UpdateGenome {
 public:
  explicit UpdateGenome(const Params& params)
      : n_loc_(params.geno.n_loc),
        n_ind_(params.geno.n_ind),
        n_sex_(n_ind_ / 2),
        v_rec_(params.geno.v_rec),
        v_mut_(params.geno.v_mut),
        bw_() {};
  void operator()(State& state);

 private:
  const std::size_t n_loc_;
  const std::size_t n_ind_;
  const std::size_t n_sex_;
  const Eigen::VectorXd& v_rec_;
  const Eigen::VectorXd& v_mut_;
  const double* ptr_rec_;
  const double* ptr_mut_;
  std::vector<std::uint64_t> transmit_chunk_;
  rng::BernoulliWord<16> bw_;

  std::uint64_t gamWord(
      std::uint64_t h0, std::uint64_t h1, std::size_t valid = 64);

  void updateGenome(State& state);
  void updateNurture(State& state);
};

inline std::uint64_t UpdateGenome::gamWord(
    std::uint64_t ind_h0, std::uint64_t ind_h1, std::size_t valid) {
  // set recombination probabilities for loci in word
  bw_.set_probs(ptr_rec_, valid);

  // sample a 0-1 recombination mask
  std::uint64_t par = bw_.sample();

  // select the initial parental strand uniformly
  bool par0 = bw_.coinflip();

  // Hallis-Steele shift cumulative sum mod 2
  par ^= par << 1;
  par ^= par << 2;
  par ^= par << 4;
  par ^= par << 8;
  par ^= par << 16;
  par ^= par << 32;
  if (par0) par = ~par;

  // set mutation probabilities for the loci
  bw_.set_probs(ptr_mut_, valid);
  std::uint64_t mut = bw_.sample();

  return ((par & ind_h0) | (~par & ind_h1)) ^ mut;
}

inline void UpdateGenome::operator()(State& state) {
  if (state.geno().view() != HaploView::IndividualMajor)
    throw std::runtime_error("update genome requires ind-major view");

  constexpr std::size_t IncWord = 64;
  const std::size_t n_words = state.geno().n_words();

  HaploBuf& h0 = state.geno().h0();
  HaploBuf& h1 = state.geno().h1();
  HaploBuf& h0_off = state.geno(Generation::Parents).h0();
  HaploBuf& h1_off = state.geno(Generation::Parents).h1();
  const Matching& matching = state.matching();

  for (std::size_t pair = 0; pair < n_sex_; ++pair) {
    std::size_t fpair = matching[pair] + n_sex_;
    ptr_rec_ = v_rec_.data();
    ptr_mut_ = v_mut_.data();
    std::size_t valid = 64;

    for (std::size_t word = 0; word < n_words; ++word) {
      if (word == n_words - 1)
        valid = (n_loc_ % 64 == 0) ? 64 : n_loc_ % 64;

      std::uint64_t male_h0 = h0(pair, word);
      std::uint64_t male_h1 = h1(pair, word);
      std::uint64_t female_h0 = h0(fpair, word);
      std::uint64_t female_h1 = h1(fpair, word);

      // male child
      h0_off(pair, word) = gamWord(male_h0, male_h1, valid);
      h1_off(pair, word) = gamWord(female_h0, female_h1, valid);

      // female child
      h0_off(fpair, word) = gamWord(male_h0, male_h1, valid);
      h1_off(fpair, word) = gamWord(female_h0, female_h1, valid);

      ptr_rec_ += IncWord;
      ptr_mut_ += IncWord;
    }
  }
}

}  // namespace amsim
