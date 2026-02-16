#pragma once

#include <amsim/rng.h>
#include <amsim/utils.h>
#include <amsim/params.h>

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace amsim::genome {

enum class HaploView : bool { LocusMajor, IndividualMajor };

class HaploBuf {
 public:
  explicit HaploBuf(std::size_t n_ind, std::size_t n_loc);

  std::size_t n_ind() const noexcept { return n_ind_; }
  std::size_t n_loc() const noexcept { return n_loc_; }
  std::size_t n_rows() const noexcept { return n_rows_; }
  std::size_t n_words() const noexcept { return n_words_; }

  HaploView view() const noexcept { return view_; }

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
  std::size_t n_words_;
  std::vector<std::uint64_t> buf_;
  HaploView view_;
};

class GenoBuf {
 public:
  explicit GenoBuf(const Params& params);

  Eigen::VectorXd& v_lmean() noexcept { return v_lmean_; }
  Eigen::VectorXd& v_lvar() noexcept { return v_lvar_; }
  Eigen::VectorXd& v_lmaf() noexcept { return v_lmaf_; }

  double v_lmean(std::size_t loc) const noexcept { return v_lmean_(loc); }
  double v_lvar(std::size_t loc) const noexcept { return v_lvar_(loc); }
  double v_lmaf(std::size_t loc) const noexcept { return v_lmaf_(loc); }

  std::size_t n_ind() const noexcept { return h0_.n_ind(); }
  std::size_t n_loc() const noexcept { return h0_.n_loc(); }
  std::size_t n_rows() const noexcept { return h0_.n_rows(); }
  std::size_t n_words() const noexcept { return h0_.n_words(); }

  HaploBuf& h0() noexcept { return h0_; }
  HaploBuf& h1() noexcept { return h1_; }
  const HaploBuf& h0() const noexcept { return h0_; }
  const HaploBuf& h1() const noexcept { return h1_; }
  HaploView view() const noexcept { return h0_.view(); }

  std::array<std::uint64_t, 2> gam_word(
      std::uint64_t ind_h0,
      std::uint64_t ind_h1,
      const double* v_rec,
      const double* v_mut) noexcept;

  void transpose() noexcept;
  void compute_mafs();
  void compute_stats();

  void decompress(
      std::size_t ind_start,
      std::size_t ind_end,
      const std::vector<std::size_t>& loc,
      Eigen::MatrixXd& out,
      bool standardise = false);

 private:
  const Eigen::VectorXd v_mut_;
  const Eigen::VectorXd v_rec_;
  const Eigen::VectorXd v_maf_;
  Eigen::VectorXd v_lmean_;
  Eigen::VectorXd v_lvar_;
  Eigen::VectorXd v_lmaf_;
  rng::BernoulliWord<16> bw_;
  HaploBuf h0_;
  HaploBuf h1_;
};

}  // namespace amsim::genome
