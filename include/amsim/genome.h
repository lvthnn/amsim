#ifndef AMSIMCPP_GENOME_H
#define AMSIMCPP_GENOME_H

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

  std::vector<double>& v_lmean() noexcept { return v_lmean_; }
  std::vector<double>& v_lvar() noexcept { return v_lvar_; }
  std::vector<double>& v_lmaf() noexcept { return v_lmaf_; }
  double v_lmean(std::size_t loc) const noexcept { return v_lmean_[loc]; }
  double v_lvar(std::size_t loc) const noexcept { return v_lvar_[loc]; }
  double v_lmaf(std::size_t loc) const noexcept { return v_lmaf_[loc]; }
  std::size_t n_loc() const noexcept { return h0_.n_loc(); }
  std::size_t n_ind() const noexcept { return h0_.n_ind(); }
  HaploView view() const noexcept { return h0_.view(); }
  HaploBuf& h0() noexcept { return h0_; }
  HaploBuf& h1() noexcept { return h1_; }
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
  std::vector<double> v_lvar_;
  std::vector<double> v_lmaf_;
  rng::BW16 bw_;
  HaploBuf h0_;
  HaploBuf h1_;
  uint64_t gamWord(std::uint64_t ind_h0, std::uint64_t ind_h1) noexcept;
};

}  // namespace amsim

#endif  // AMSIMCPP_GENOME_H
