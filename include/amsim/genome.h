#ifndef AMSIMCPP_GENOME_H
#define AMSIMCPP_GENOME_H

#pragma once

#include <amsim/haplobuf.h>
#include <amsim/haploview.h>
#include <amsim/rng.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace amsim {

class Genome {
 public:
  Genome(
      std::size_t n_ind,
      std::size_t n_loc,
      std::vector<double> v_mut,
      std::vector<double> v_rec,
      std::vector<double> v_maf,
      const rng::Xoshiro256ss& rng);

  inline std::vector<double>& v_lmean() noexcept { return v_lmean_; }
  inline std::vector<double>& v_lvar() noexcept { return v_lvar_; }
  inline std::vector<double>& v_lmaf() noexcept { return v_lmaf_; }
  inline double v_lmean(std::size_t loc) const noexcept {
    return v_lmean_[loc];
  }
  inline double v_lvar(std::size_t loc) const noexcept { return v_lvar_[loc]; }
  inline double v_lmaf(std::size_t loc) const noexcept { return v_lmaf_[loc]; }
  inline std::size_t n_loc() const noexcept { return H0_.n_loc(); }
  inline std::size_t n_ind() const noexcept { return H0_.n_ind(); }
  inline HaploView view() const noexcept { return H0_.view(); }
  inline HaploBuf& H0() noexcept { return H0_; }
  inline HaploBuf& H1() noexcept { return H1_; }
  void generate_haplotypes() noexcept;
  void transpose() noexcept;
  void compute_mafs();
  void compute_stats();
  void update(std::vector<std::size_t> matching);

 private:
  const std::vector<double> v_mut_;
  const std::vector<double> v_rec_;
  const std::vector<double> v_maf_;
  std::vector<double> v_lmean_;
  std::vector<double> v_bmean;
  std::vector<double> v_lvar_;
  std::vector<double> v_lmaf_;
  rng::BW16 bw_;
  HaploBuf H0_;
  HaploBuf H1_;
  uint64_t gam_word_(std::uint64_t ind_H0, std::uint64_t ind_H1) noexcept;
};

}  // namespace amsim

#endif  // AMSIMCPP_GENOME_H
