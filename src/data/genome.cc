#include <amsim/data/genome.h>
#include <amsim/params.h>
#include <amsim/utils.h>

#include <cstddef>
#include <cstdint>

namespace amsim::genome {

HaploBuf::HaploBuf(std::size_t n_ind, std::size_t n_loc)
    : n_ind_(n_ind),
      n_loc_(n_loc),
      n_rows_((n_loc_ + 63) & ~static_cast<std::size_t>(63)),
      n_words_((n_ind_ + 63) / 64),
      view_(HaploView::LocusMajor) {
  buf_.resize(n_rows_ * n_words_);
}

void HaploBuf::transpose() noexcept {
  // pre-transpose
  const std::size_t src_rows = n_rows_;
  const std::size_t src_words = n_words_;

  // post-transpose
  const std::size_t dst_rows =
      (view_ == HaploView::LocusMajor)
          ? ((n_ind_ + 63) & ~static_cast<std::size_t>(63))
          : ((n_loc_ + 63) & ~static_cast<std::size_t>(63));

  const std::size_t dst_cols_w = (view_ == HaploView::LocusMajor)
                                     ? ((n_loc_ + 63) / 64)
                                     : ((n_ind_ + 63) / 64);

  // allocate destination buffer
  std::vector<std::uint64_t> out(dst_rows * dst_cols_w, 0);

  // src_rows is multiple of 64 by construction
  const std::size_t src_n_rows = src_rows / 64;
  const std::size_t src_n_cols = src_words;

  std::uint64_t tile[64];

  for (std::size_t ct = 0; ct < src_n_rows; ++ct) {
    for (std::size_t cb = 0; cb < src_n_cols; ++cb) {
      for (std::size_t r = 0; r < 64; ++r) tile[r] = (*this)((ct * 64) + r, cb);

      utils::bitmatrix_transpose(tile);

      for (std::size_t r = 0; r < 64; ++r) {
        const std::size_t dst_row = (cb * 64) + r;
        const std::size_t dst_col = ct;
        out[(dst_row * dst_cols_w) + dst_col] = tile[r];
      }
    }
  }

  buf_.swap(out);
  n_rows_ = dst_rows;
  n_words_ = dst_cols_w;
  view_ = (view_ == HaploView::LocusMajor) ? HaploView::IndividualMajor
                                           : HaploView::LocusMajor;
}

GenoBuf::GenoBuf(const Params& params)
    : v_mut_(std::move(params.geno.v_mut)),
      v_rec_(std::move(params.geno.v_rec)),
      v_maf_(std::move(params.geno.v_maf)),
      v_lmean_(params.geno.n_loc),
      v_lvar_(params.geno.n_loc),
      v_lmaf_(params.geno.n_loc),
      bw_(),
      h0_(params.geno.n_ind, params.geno.n_loc),
      h1_(params.geno.n_ind, params.geno.n_loc) {};

std::array<std::uint64_t, 2> GenoBuf::gam_word(
    std::uint64_t ind_h0,
    std::uint64_t ind_h1,
    const double* v_rec,
    const double* v_mut) noexcept {
  // Set recombination probabilities for loci in word
  bw_.set_probs(v_rec);

  // Sample a 0-1 recombination mask
  std::uint64_t par = bw_.sample();

  // Select the initial parental strand uniformly
  bool par0 = bw_.coinflip();

  // Hallis-Steele shift cumulative sum mod 2
  par ^= par << 1;
  par ^= par << 2;
  par ^= par << 4;
  par ^= par << 8;
  par ^= par << 16;
  par ^= par << 32;
  if (par0) par = ~par;

  // Set mutation probabilities for the loci
  bw_.set_probs(v_mut);
  std::uint64_t mut = bw_.sample();

  return {par, ((par & ind_h0) | (~par & ind_h1)) ^ mut};
}

void GenoBuf::transpose() noexcept {
  h0_.transpose();
  h1_.transpose();
}

void GenoBuf::compute_mafs() {
  if (h0_.view() != HaploView::LocusMajor)
    throw std::runtime_error(
        "GenoBuf::compute_mafs: compute MAFs in locus-major view.");

  std::size_t n_ind = h0_.n_ind();
  std::size_t n_loc = h0_.n_loc();

  const std::size_t n_bloc_ind = (n_ind + 63) / 64;

  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    // Total locus dosage
    std::size_t ct_loc = 0;

    for (std::size_t bloc = 0; bloc < n_bloc_ind; ++bloc) {
      if (bloc == n_bloc_ind - 1 && (n_ind % 64)) {
        // Mask to trim off padding in last block
        std::uint64_t mask = (1ULL << (n_ind % 64)) - 1ULL;
        ct_loc += __builtin_popcountll(h0_(loc, bloc) & mask);
        ct_loc += __builtin_popcountll(h1_(loc, bloc) & mask);
      } else {
        // Popcount words for fast total haplotype dosage
        ct_loc += __builtin_popcountll(h0_(loc, bloc));
        ct_loc += __builtin_popcountll(h1_(loc, bloc));
      }
    }

    // Compute the MAF
    v_lmaf_(loc) =
        static_cast<double>(ct_loc) / (2.0 * static_cast<double>(n_ind));
  }
}

void GenoBuf::compute_stats() {
  if (h0_.view() != HaploView::LocusMajor)
    throw std::runtime_error(
        "GenoBuf::compute_stats: compute stats in loc-major view.");

  std::size_t n_ind = h0_.n_ind();
  std::size_t n_loc = h0_.n_loc();

  const std::size_t n_bloc_ind = (n_ind + 63) / 64;

  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    std::size_t hom = 0;
    std::size_t het = 0;
    for (std::size_t bloc = 0; bloc < n_bloc_ind; ++bloc) {
      if (bloc == n_bloc_ind - 1 && (n_ind % 64)) {
        std::uint64_t mask = (1ULL << (n_ind % 64)) - 1ULL;
        hom += __builtin_popcountll((h0_(loc, bloc) & h1_(loc, bloc)) & mask);
        het += __builtin_popcountll((h0_(loc, bloc) ^ h1_(loc, bloc)) & mask);
      } else {
        hom += __builtin_popcountll(h0_(loc, bloc) & h1_(loc, bloc));
        het += __builtin_popcountll(h0_(loc, bloc) ^ h1_(loc, bloc));
      }
    }
    double mean_loc = static_cast<double>((2.0 * hom) + het) / n_ind;
    double mean_sqloc = static_cast<double>((4.0 * hom) + het) / n_ind;
    double var_loc = mean_sqloc - (mean_loc * mean_loc);

    v_lmean_(loc) = mean_loc;
    v_lvar_(loc) = var_loc;
  }
}

void GenoBuf::decompress(
    std::size_t ind_start,
    std::size_t ind_end,
    const std::vector<std::size_t>& loc,
    Eigen::MatrixXd& out,
    bool standardise) {
  if (view() != HaploView::LocusMajor)
    throw std::runtime_error("GenoBuf::decompress: require loc-major view");
  if (ind_end > n_ind())
    throw std::runtime_error(
        "GenoBuf::decompress: argument ind_end exceeds number of individuals");

  std::size_t n_loc = loc.size();
  std::size_t start = ind_start / 64;
  std::size_t end = (ind_end + 63) / 64;

  for (std::size_t el = 0; el < n_loc; ++el) {
    double scl = 1.0 / std::sqrt(v_lvar_[loc[el]]);
    double cen = -v_lmean_[loc[el]] * scl;
    for (std::size_t word = start; word < end; ++word) {
      std::size_t bit_lo = (word == start) ? (ind_start % 64) : 0;
      std::size_t bit_hi = (word == end - 1) ? ((ind_end - 1) % 64) + 1 : 64;
      std::uint64_t h0 = h0_(loc[el], word);
      std::uint64_t h1 = h1_(loc[el], word);

      for (std::size_t bit = bit_lo; bit < bit_hi; ++bit) {
        out((64 * word) + bit - ind_start, el) =
            (standardise)
                ? (scl * (((h0 >> bit) & 1ULL) + ((h1 >> bit) & 1ULL))) + cen
                : ((h0 >> bit) & 1ULL) + ((h1 >> bit) & 1ULL);
      }
    }
  }
}

}  // namespace amsim::genome
