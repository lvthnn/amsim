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

#include <amsim/core/params.h>
#include <amsim/core/rng.h>
#include <amsim/core/utils.h>

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace amsim {

inline void bitmatrixSwap(
    std::uint64_t matrix[], std::size_t width, std::uint64_t mask) {
  std::size_t inner;
  std::size_t outer;
  for (outer = 0; outer < 64 / (width * 2); ++outer) {
    for (inner = 0; inner < width; ++inner) {
      std::uint64_t* x = &matrix[(inner) + (outer * width * 2)];
      std::uint64_t* y = &matrix[(inner + width) + (outer * width * 2)];
      *x = ((*y << width) & mask) ^ *x;
      *y = ((*x & mask) >> width) ^ *y;
      *x = ((*y << width) & mask) ^ *x;
    }
  }
}

inline void bitmatrixTranspose(std::uint64_t* matrix) {
  std::size_t swap_width = 64;
  auto swap_mask = static_cast<std::uint64_t>(-1);
  while (swap_width != 1) {
    swap_width >>= 1;
    swap_mask = swap_mask ^ (swap_mask >> swap_width);
    bitmatrixSwap(matrix, swap_width, swap_mask);
  }
}

enum class BufferLayout : bool { LocusMajor, IndividualMajor };

class HaplotypeBuffer {
 public:
  HaplotypeBuffer(std::size_t n_individuals, std::size_t n_loc)
      : n_ind_(n_individuals),
        n_loc_(n_loc),
        n_rows_((n_loc_ + 63) & ~static_cast<std::size_t>(63)),
        n_words_((n_ind_ + 63) / 64),
        buf_(n_rows_ * n_words_) {}

  std::size_t numIndividuals() const noexcept { return n_ind_; }
  std::size_t numLoci() const noexcept { return n_loc_; }
  std::size_t numRows() const noexcept { return n_rows_; }
  std::size_t numWords() const noexcept { return n_words_; }

  BufferLayout layout() const noexcept { return layout_; }

  std::uint64_t& operator()(std::size_t i, std::size_t j) noexcept {
    return buf_[(i * n_words_) + j];
  }

  const std::uint64_t& operator()(std::size_t i, std::size_t j) const noexcept {
    return buf_[(i * n_words_) + j];
  }

  std::uint64_t* rowptr(std::size_t i) noexcept { return &buf_[i * n_words_]; }

  const std::uint64_t* rowptr(std::size_t i) const noexcept {
    return &buf_[i * n_words_];
  }

  void transpose() noexcept;

 private:
  const std::size_t n_ind_;
  const std::size_t n_loc_;
  std::size_t n_rows_;
  std::size_t n_words_{};
  std::vector<std::uint64_t> buf_;
  BufferLayout layout_{};
};

inline void HaplotypeBuffer::transpose() noexcept {
  // pre-transpose
  const std::size_t src_rows = n_rows_;
  const std::size_t src_words = n_words_;

  // post-transpose
  const std::size_t dst_rows =
      (layout_ == BufferLayout::LocusMajor)
          ? ((n_ind_ + 63) & ~static_cast<std::size_t>(63))
          : ((n_loc_ + 63) & ~static_cast<std::size_t>(63));

  const std::size_t dst_cols_w = (layout_ == BufferLayout::LocusMajor)
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

      bitmatrixTranspose(tile);

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
  layout_ = (layout_ == BufferLayout::LocusMajor)
                ? BufferLayout::IndividualMajor
                : BufferLayout::LocusMajor;
}

class GenotypeBuffer {
 public:
  explicit GenotypeBuffer(const Params& params)
      : loc_mut_(std::move(params.geno.locus_mut)),
        loc_rec_(std::move(params.geno.locus_rec)),
        loc_freq_(std::move(params.geno.locus_freq)),
        loc_lmean_(params.geno.n_loc),
        loc_lvar_(params.geno.n_loc),
        loc_lfreq_(params.geno.n_loc),
        bw_(),
        h0_(params.global.n_ind, params.geno.n_loc),
        h1_(params.global.n_ind, params.geno.n_loc) {};

  Eigen::VectorXd& locusMean() noexcept { return loc_lmean_; }
  Eigen::VectorXd& locusVar() noexcept { return loc_lvar_; }
  Eigen::VectorXd& locusFreq() noexcept { return loc_lfreq_; }

  const Eigen::VectorXd& locusMean() const noexcept { return loc_lmean_; }
  const Eigen::VectorXd& locusVar() const noexcept { return loc_lvar_; }
  const Eigen::VectorXd& locusFreq() const noexcept { return loc_lfreq_; }

  double locusMean(std::size_t loc) const noexcept { return loc_lmean_(loc); }
  double locusVar(std::size_t loc) const noexcept { return loc_lvar_(loc); }
  double locusFreq(std::size_t loc) const noexcept { return loc_lfreq_(loc); }

  std::size_t numIndividuals() const noexcept { return h0_.numIndividuals(); }
  std::size_t numLoci() const noexcept { return h0_.numLoci(); }
  std::size_t numRows() const noexcept { return h0_.numRows(); }
  std::size_t numWords() const noexcept { return h0_.numWords(); }

  HaplotypeBuffer& h0() noexcept { return h0_; }
  HaplotypeBuffer& h1() noexcept { return h1_; }
  const HaplotypeBuffer& h0() const noexcept { return h0_; }
  const HaplotypeBuffer& h1() const noexcept { return h1_; }
  BufferLayout view() const noexcept { return h0_.layout(); }

  void transpose() noexcept;
  void computeLocusFreqs();
  void computeLocusStats();

  void decompress(
      std::size_t ind_start,
      std::size_t ind_end,
      const std::vector<std::size_t>& loc,
      Eigen::MatrixXd& out,
      bool centre = false,
      bool scale = false);

 private:
  const Eigen::VectorXd loc_mut_;
  const Eigen::VectorXd loc_rec_;
  const Eigen::VectorXd loc_freq_;
  Eigen::VectorXd loc_lmean_;
  Eigen::VectorXd loc_lvar_;
  Eigen::VectorXd loc_lfreq_;
  rng::BernoulliWord<16> bw_;
  HaplotypeBuffer h0_;
  HaplotypeBuffer h1_;
};

inline void GenotypeBuffer::transpose() noexcept {
  h0_.transpose();
  h1_.transpose();
}

inline void GenotypeBuffer::computeLocusFreqs() {
  if (h0_.layout() != BufferLayout::LocusMajor)
    throw std::runtime_error(
        "GenoBuf::compute_mafs: compute MAFs in locus-major view.");

  std::size_t n_ind = h0_.numIndividuals();
  std::size_t n_loc = h0_.numLoci();

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
    loc_lfreq_(loc) =
        static_cast<double>(ct_loc) / (2.0 * static_cast<double>(n_ind));
  }
}

inline void GenotypeBuffer::computeLocusStats() {
  if (h0_.layout() != BufferLayout::LocusMajor)
    throw std::runtime_error(
        "GenoBuf::compute_stats: compute stats in loc-major view.");

  std::size_t n_ind = h0_.numIndividuals();
  std::size_t n_loc = h0_.numLoci();

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

    loc_lmean_(loc) = mean_loc;
    loc_lvar_(loc) = var_loc;
  }
}

inline void GenotypeBuffer::decompress(
    std::size_t ind_start,
    std::size_t ind_end,
    const std::vector<std::size_t>& loc,
    Eigen::MatrixXd& out,
    bool centre,
    bool scale) {
  if (view() != BufferLayout::LocusMajor)
    throw std::runtime_error("GenoBuf::decompress: require loc-major view");
  if (ind_end > numIndividuals())
    throw std::runtime_error(
        "GenBuf::decompress: argument ind_end exceeds number of individuals");

  std::size_t n_loc = loc.size();
  std::size_t start = ind_start / 64;
  std::size_t end = (ind_end + 63) / 64;

  for (std::size_t el = 0; el < n_loc; ++el) {
    double scl = (scale) ? 1.0 / std::sqrt(loc_lvar_[loc[el]]) : 1.0;
    double cen = (centre) ? -loc_lmean_[loc[el]] * scl : 0;

    for (std::size_t word = start; word < end; ++word) {
      std::size_t bit_lo = (word == start) ? (ind_start % 64) : 0;
      std::size_t bit_hi = (word == end - 1) ? ((ind_end - 1) % 64) + 1 : 64;
      std::uint64_t h0 = h0_(loc[el], word);
      std::uint64_t h1 = h1_(loc[el], word);

      for (std::size_t bit = bit_lo; bit < bit_hi; ++bit) {
        out((64 * word) + bit - ind_start, el) =
            (scl * (((h0 >> bit) & 1ULL) + ((h1 >> bit) & 1ULL))) + cen;
      }
    }
  }
}

}  // namespace amsim
