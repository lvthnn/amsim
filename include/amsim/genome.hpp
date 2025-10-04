#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include <amsim/rng.hpp>

namespace amsim {
  enum class HapMatView : bool {
    IND_MAJOR,
    LOC_MAJOR
  };
  
  class HapMat {
  public:
    HapMat(size_t n_ind, size_t n_loc);

    inline std::size_t n_ind() const noexcept { return n_ind_; }
    inline std::size_t n_loc() const noexcept { return n_loc_; }
    inline std::size_t n_rows() const noexcept { return n_rows_; }
    inline std::size_t n_words() const noexcept { return n_words_; }
    inline HapMatView view() const noexcept { return view_; }

    inline std::uint64_t& operator()(std::size_t i, std::size_t j) noexcept {
			return (view_ == HapMatView::IND_MAJOR)
				? data_[i * n_words_ + j] : data_[j * n_rows_ + i];
		};

    inline const std::uint64_t& operator()(std::size_t i, std::size_t j) const noexcept {
			return (view_ == HapMatView::IND_MAJOR)
				? data_[i * n_words_ + j] : data_[j * n_rows_ + i];
		};

    inline const std::uint64_t* rowptr(std::size_t i) const noexcept {
      if (view_ == HapMatView::IND_MAJOR) return &data_[i * n_words_];
      return nullptr;
    }

    inline std::uint64_t* rowptr(std::size_t i) noexcept {
      if (view_ == HapMatView::IND_MAJOR) return &data_[i * n_words_];
      return nullptr;
    }

    void transpose() noexcept;

  private:
    const std::size_t n_ind_;
    const std::size_t n_loc_;
    std::size_t n_rows_;
    std::size_t n_words_;
    std::vector<uint64_t> data_;
    HapMatView view_;
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
    inline HapMatView view() const noexcept { return H0_.view(); }
    inline HapMat& H0() noexcept { return H0_; }
    inline HapMat& H1() noexcept { return H1_; }
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
    HapMat H0_;
    HapMat H1_;
    uint64_t gam_word_(std::size_t ind, std::size_t word) noexcept;
  };
}
