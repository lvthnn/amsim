#include <amsim/core/params.h>
#include <amsim/core/state.h>
#include <amsim/transform.h>

namespace amsim {

std::array<std::uint64_t, 2> UpdateGenome::gamWord(
    std::uint64_t ind_h0, std::uint64_t ind_h1) {
  // set recombination probabilities for loci in word
  bw_.set_probs(ptr_rec_);

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
  bw_.set_probs(ptr_mut_);
  std::uint64_t mut = bw_.sample();

  return {par, ((par & ind_h0) | (~par & ind_h1)) ^ mut};
}

void UpdateGenome::operator()(State& state) {
  if (state.geno().view() != HaploView::IndividualMajor)
    throw std::runtime_error("update genome requires ind-major view");

  constexpr std::size_t IncWord = 64;
  const std::size_t n_words = state.geno().n_words();

  HaploBuf& h0 = state.geno().h0();
  HaploBuf& h1 = state.geno().h1();
  HaploBuf& h0_off = state.geno_par().h0();
  HaploBuf& h1_off = state.geno_par().h1();
  const Matching& matching = state.matching();

  for (std::size_t pair = 0; pair < n_sex_; ++pair) {
    std::size_t fpair = matching[pair] + n_sex_;
    ptr_rec_ = v_rec_.data();
    ptr_mut_ = v_mut_.data();
    for (std::size_t word = 0; word < n_words; ++word) {
      std::uint64_t male_h0 = h0(pair, word);
      std::uint64_t male_h1 = h1(pair, word);
      std::uint64_t female_h0 = h0(fpair, word);
      std::uint64_t female_h1 = h1(fpair, word);

      // male child
      auto [mask_mm, gam_mm] = gamWord(male_h0, male_h1);
      auto [mask_fm, gam_fm] = gamWord(female_h0, female_h1);
      h0_off(pair, word) = gam_mm;
      h1_off(pair, word) = gam_fm;

      // female child
      auto [mask_mf, gam_mf] = gamWord(male_h0, male_h1);
      auto [mask_ff, gam_ff] = gamWord(female_h0, female_h1);
      h0_off(fpair, word) = gam_mf;
      h1_off(fpair, word) = gam_ff;

      ptr_rec_ += IncWord;
      ptr_mut_ += IncWord;
    }
  }
}

}  // namespace amsim
