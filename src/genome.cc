#include <amsim/genome.h>
#include <amsim/haploview.h>
#include <amsim/logger.h>
#include <amsim/rng.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace amsim {
Genome::Genome(
    std::size_t n_ind,
    std::size_t n_loc,
    std::vector<double> v_mut,
    std::vector<double> v_rec,
    std::vector<double> v_maf,
    const rng::Xoshiro256ss& rng)
    : v_mut_(std::move(v_mut)),
      v_rec_(std::move(v_rec)),
      v_maf_(std::move(v_maf)),
      v_lmean_(n_loc),
      v_lvar_(n_loc),
      v_lmaf_(n_loc),
      bw_(rng),
      h0_(n_ind, n_loc),
      h1_(n_ind, n_loc) {};

std::uint64_t Genome::gamWord(
    std::uint64_t ind_h0, std::uint64_t ind_h1) noexcept {
  std::uint64_t par = bw_.sample();
  bool par0 = bw_.coinflip();
  if (v_rec_[0] < 0.5) {
    par ^= par << 1;
    par ^= par << 2;
    par ^= par << 4;
    par ^= par << 8;
    par ^= par << 16;
    par ^= par << 32;
    if (par0) par = ~par;
  }
  return (par & ind_h0) | (~par & ind_h1);
}

void Genome::generate_haplotypes() noexcept {
  std::size_t n_loc = h0_.n_loc();
  std::size_t n_words = h0_.n_words();

  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    bw_.set_prob(v_maf_[loc]);
    std::uint64_t* word0 = h0_.rowptr(loc);
    std::uint64_t* word1 = h1_.rowptr(loc);

    for (std::size_t word = 0; word < n_words; ++word) {
      word0[word] = bw_.sample();
      word1[word] = bw_.sample();
    }
  }
}

void Genome::transpose() noexcept {
  h0_.transpose();
  h1_.transpose();
}

void Genome::compute_mafs() {
  if (h0_.view() == HaploView::IND_MAJOR)
    throw std::runtime_error("Compute MAFs in locus-major view.");

  std::size_t n_ind = h0_.n_ind();
  std::size_t n_loc = h0_.n_loc();

  const std::size_t n_bloc_ind = (n_ind + 63) / 64;

  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    std::size_t ct_loc = 0;
    for (std::size_t bloc = 0; bloc < n_bloc_ind; ++bloc) {
      if (bloc == n_bloc_ind - 1 && (n_ind % 64)) {
        std::uint64_t mask = (1ULL << (n_ind % 64)) - 1ULL;
        ct_loc += __builtin_popcountll(h0_(loc, bloc) & mask);
        ct_loc += __builtin_popcountll(h1_(loc, bloc) & mask);
      } else {
        ct_loc += __builtin_popcountll(h0_(loc, bloc));
        ct_loc += __builtin_popcountll(h1_(loc, bloc));
      }
    }
    v_lmaf_[loc] =
        static_cast<double>(ct_loc) / (2.0 * static_cast<double>(n_ind));
  }
}

void Genome::compute_stats() {
  if (h0_.view() == HaploView::IND_MAJOR)
    throw std::runtime_error("Compute stats in loc-major view.");

  std::size_t n_ind = h0_.n_ind();
  std::size_t n_loc = h0_.n_loc();

  const std::size_t n_bloc_ind = (n_ind + 63) / 64;

  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    std::size_t hom = 0;
    std::size_t het = 0;
    for (std::size_t bloc = 0; bloc < n_bloc_ind; ++bloc) {
      hom += __builtin_popcountll(h0_(loc, bloc) & h1_(loc, bloc));
      het += __builtin_popcountll(h0_(loc, bloc) ^ h1_(loc, bloc));
    }
    double mean_loc = static_cast<double>((2.0 * hom) + het) / n_ind;
    double mean_sqloc = static_cast<double>((4.0 * hom) + het) / n_ind;
    double var_loc = mean_sqloc - (mean_loc * mean_loc);

    v_lmean_[loc] = mean_loc;
    v_lvar_[loc] = var_loc;
  }
}

void Genome::update(std::vector<std::size_t> matching) {
  if (h0_.view() == HaploView::LOC_MAJOR)
    throw std::runtime_error("Update in ind-major view.");

  const std::size_t n_ind = h0_.n_ind();
  const std::size_t n_words = h0_.n_words();
  const std::size_t n_pairs = n_ind / 2;
  bw_.set_prob(v_rec_[0]);

  // @TODO: integrate mutation
  for (std::size_t pair = 0; pair < n_pairs; ++pair) {
    std::size_t fpair = matching[pair] + n_pairs;
    for (std::size_t word = 0; word < n_words; ++word) {
      std::uint64_t male_h0 = h0_(pair, word);
      std::uint64_t male_h1 = h1_(pair, word);
      std::uint64_t female_h0 = h0_(fpair, word);
      std::uint64_t female_h1 = h1_(fpair, word);

      // male child
      h0_(pair, word) = gamWord(male_h0, male_h1);
      h1_(pair, word) = gamWord(female_h0, female_h1);

      // female child
      h0_(fpair, word) = gamWord(male_h0, male_h1);
      h1_(fpair, word) = gamWord(female_h0, female_h1);
    }
  }
}
}  // namespace amsim
