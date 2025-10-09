//------------------------------------------------------------------------------
// amsimcpp : phenoarch.h
//------------------------------------------------------------------------------

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <vector>
#include <random>
#include <numeric>
#include <cassert>
#include <iostream>

#include <amsim/phenoarch.h>
#include <amsim/rng.h>

#if defined(__APPLE__)
  #include<Accelerate/Accelerate.h>
#elif defined(USE_BLAS)
  #include<cblas.h>
  #include<lapacke.h>
#endif

namespace amsim {
  PhenoArch::PhenoArch(std::size_t n_pheno, std::size_t n_loc_total,
                       std::vector<std::size_t> n_loc, std::vector<double> h2_gen,
                       std::vector<double> gen_cor, std::vector<double> env_cor,
                       const rng::Xoshiro256ss &rng)
    : n_pheno_(n_pheno),
      n_loc_tot_(n_loc_total),
      n_loc_(std::move(n_loc)),
      h2_gen_(std::move(h2_gen)),
      gen_cor_(std::move(gen_cor)),
      env_chol_(std::move(env_cor)),
      rng_polar_(rng),
      rng_unf_(rng)  {
    assert(env_cor.size() == n_pheno_ * n_pheno_);
    assert(h2_gen_.size() == n_pheno_ * n_pheno_);
    assert(n_loc_.size() == n_pheno_ * n_pheno_);
    assert(gen_cor_.size() == n_pheno_ * n_pheno_);

    env_chol_.resize(n_pheno_ * n_pheno_);

    // cast to LAPACK-legible form
    char clpk_uplo_ = 'L';
    int clpk_n_pheno_ = static_cast<int>(n_pheno_);
    int clpk_lda_ = 3;
    int clpk_out_;

    #if defined(__APPLE__) && defined(USE_BLAS)
      dpotrf_(&clpk_uplo_, &clpk_n_pheno_, env_chol_.data(), &clpk_lda_,
              &clpk_out_);
    #elif defined(__linux__) && defined(USE_BLAS)
      LAPACKE_dpotrf(&clpk_uplo_, &clpk_n_pheno_, env_chol_.data(), &clpk_lda_,
                     &clpk_out_);
    #endif

    for (std::size_t c = 0; c < n_pheno; c++)
      for (std::size_t r = 0; r < c; r++)
        env_chol_[c * n_pheno_ + r] = 0.0;
  }

  void PhenoArch::gen_env(double *ptr_env, const std::size_t n_ind) {
    rng_polar_.fill(ptr_env, n_ind * n_pheno_);
    #if defined(USE_BLAS)
      cblas_dtrmm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                  CblasNonUnit, n_pheno_, n_ind, 1, env_chol_.data(), n_pheno_,
                  ptr_env, n_pheno_);
    #else
      throw std::runtime_error("Not implemented");
    #endif
  }

  std::vector<double> PhenoArch::init_weights_() const {
    std::vector<double> cost_vec((n_pheno_ * (n_pheno_ + 1)) / 2);
    std::size_t idx = 0;
    for (std::size_t i = 0; i < n_pheno_; i++) {
      for (std::size_t j = i; j < n_pheno_; j++) {
        cost_vec[idx] = gen_cor_[j * n_pheno_ + i] * std::sqrt(n_loc_[i] * n_loc_[j]);
        idx++;
      }
    }
    return cost_vec;
  }

  std::vector<uint64_t> PhenoArch::init_mask_() {
    const std::size_t n_words = (n_loc_tot_ + 63) / 64;
    std::vector<uint64_t> loc_mask(n_words * n_pheno_);

    for (std::size_t pheno = 0; pheno < n_pheno_; pheno++) {
      uint64_t* loc_ptr = &loc_mask[pheno * n_words];
      const std::size_t n_loc_pheno = n_loc_[pheno];
      for (std::size_t r_id = n_loc_tot_ - n_loc_pheno; r_id < n_loc_tot_; r_id++) {
        const std::size_t l_id = rng_unf_.sample(r_id + 1);

        // @TODO: Change this so we write in row-major format (makes optimisation step less computationally heavy)
        const std::size_t lw = l_id / 64; const std::size_t lo = l_id % 64;
        const std::size_t rw = r_id / 64; const std::size_t ro = r_id % 64;

        if (loc_ptr[lw] & (1ull << lo)) loc_ptr[rw] |= (1ull << ro);
        else loc_ptr[lw] |= (1ull << lo);
      }
    }
    return loc_mask;
  }

  std::vector<std::size_t> PhenoArch::init_intersect_(const std::vector<uint64_t> &mask) const {
    std::vector<std::size_t> intersect(n_pheno_ * (n_pheno_ - 1) / 2);
    const std::size_t n_words = (n_loc_tot_ + 63) / 64;

    std::size_t id = 0;
    for (std::size_t i = 0; i < n_pheno_; i++) {
      for (std::size_t j = i + 1; j < n_pheno_; j++) {
        const uint64_t* pheno_i = &mask[n_words * i];
        const uint64_t* pheno_j = &mask[n_words * j];
        for (std::size_t word = 0; word < n_words; word++) {
          intersect[id] += __builtin_popcountll(pheno_i[word] & pheno_j[word]);
        }
        id++;
      }
    }
    return intersect;
  }

  // void PhenoArch::optim_arch(double eps, std::size_t max_it) {
  //   const std::size_t n_words = (n_loc_tot_ + 63) / 64;
  //   std::vector<double> init_weights = init_weights_();
  //   std::vector<uint64_t> mask = init_mask_();

  //   for (std::size_t it = 0; it < max_it; it++) {
  //     // index of the phenotype currently inspected
  //     std::size_t pheno = it % n_pheno_;

  //     // pointer to phenotype locus mask in buffer
  //     uint64_t *ptr_pheno = &mask[n_words * pheno];

  //     for (std::size_t loc = 0; loc < n_loc_tot_; loc++) {
  //       // block containing locus in mask
  //       std::size_t block = loc / 64;

  //       // relative position of locus in the mask
  //       std::size_t offset = loc % 64;

  //       // what kind of operation are we considering?
  //       ArchOp op = (ptr_pheno[block] & (1ull << offset))
  //         ? ArchOp::ADDITION
  //         : ArchOp::DELETION;

  //       // scan over same locus in the other phenotypes
  //       for (std::size_t pheno_adj = 0; pheno_adj < n_pheno_; pheno_adj++) {
  //         if (pheno == pheno_adj) continue;

  //         // pointer to adjacent phenotype locus mask in buffer
  //         uint64_t *ptr_pheno_adj = &mask[n_words * pheno_adj];

  //         // if the relevant site is activated, increment or decrement the
  //         // intersection list depending on
  //       }

  //     }
  //   }
  // }
}