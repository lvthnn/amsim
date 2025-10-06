#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include <amsim/rng.h>

namespace amsim {
  enum class HaploView : bool {
    IND_MAJOR,
    LOC_MAJOR
  };
  
  class HaploBuf {
  public:
    HaploBuf(size_t n_ind, size_t n_loc);

    inline std::size_t n_ind() const noexcept { return n_ind_; }
    inline std::size_t n_loc() const noexcept { return n_loc_; }
    inline std::size_t n_rows() const noexcept { return n_rows_; }
    inline std::size_t n_words() const noexcept { return n_words_; }
    inline HaploView view() const noexcept { return view_; }

    inline std::uint64_t& operator()(std::size_t i, std::size_t j) noexcept {
      return data_[i * n_words_ + j];
		}

    inline const std::uint64_t& operator()(std::size_t i, std::size_t j) const noexcept {
      return data_[i * n_words_ + j];
		}

    inline const std::uint64_t* rowptr(std::size_t i) const noexcept {
      return &data_[i * n_words_];
    }

    inline std::uint64_t* rowptr(std::size_t i) noexcept {
      return &data_[i * n_words_];
    }

    void transpose() noexcept;

  private:
    const std::size_t n_ind_;
    const std::size_t n_loc_;
    std::size_t n_rows_;
    std::size_t n_words_;
    std::vector<uint64_t> data_;
    HaploView view_;
  };

  class Genome {
  public:
    Genome(size_t n_ind, size_t n_loc, std::vector<double> v_mut,
           std::vector<double> v_rec, std::vector<double> v_maf,
           uint64_t rng_seed);
    inline std::vector<double>& v_lmean() noexcept { return v_lmean_; }
    inline std::vector<double>& v_lvar() noexcept { return v_lvar_; }
    inline std::vector<double>& v_lmaf() noexcept { return v_lmaf_; }
    inline double v_lmean(std::size_t loc) const noexcept { return v_lmean_[loc]; }
    inline double v_lvar(std::size_t loc) const noexcept { return v_lvar_[loc]; }
    inline double v_lmaf(std::size_t loc) const noexcept { return v_lmaf_[loc]; }
    inline HaploView view() const noexcept { return H0_.view(); }
    inline HaploBuf& H0() noexcept { return H0_; }
    inline HaploBuf& H1() noexcept { return H1_; }
    void generate_haplotypes() noexcept;
    void transpose();
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
    uint64_t gam_word_(std::size_t ind, std::size_t word) noexcept;
  };
}
