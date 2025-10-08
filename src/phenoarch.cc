//------------------------------------------------------------------------------
// amsimcpp : phenoarch.h
//------------------------------------------------------------------------------

#include <cstddef>
#include <cstdint>
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
  PhenoArch::PhenoArch(std::size_t n_pheno, std::vector<std::size_t> n_loc,
                       std::size_t n_loc_total, std::vector<double> h2_gen,
                       std::vector<double> gen_cor, std::vector<double> env_cor,
                       const rng::Xoshiro256ss &rng)
    : n_pheno_(n_pheno),
      n_loc_(std::move(n_loc)),
      n_loc_tot_(n_loc_total),
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

  void PhenoArch::optim_arch(double eps, std::size_t max_it) {
    /**
     * declare a bit-packed uint64_t buffer for activation masks
     *
     * while (|E[cors] - cors_input| > eps) && (it < max_it)
     *   for each column
     *     for each word in bit-packed mask
     *       for each bit
     *         compute the effect of flipping on the cost
     *         if minimal remove -> update min_remove
     *         if maximal add -> update max_add
     *     engage max_add bit and disengage min_remove bit
     *
     * convert masks to std::size_t vectors for each phenotype
     *
     * allow phenotypes to access the masks
     */
    // @TODO: Implement feasibility checking of inputs given constraints
    const std::size_t n_words = (n_loc_tot_ + 63) / 64;
    std::vector<uint64_t> loc_mask(n_pheno_ * n_words);
    std::vector<std::size_t> cnt_pheno(n_pheno_);

    // generate random-state bit-masks for each phenotype using Floyd's sampling
    // without replacement algorithm
    for (std::size_t pheno = 0; pheno < n_pheno_; pheno++) {
      uint64_t* loc_ptr = &loc_mask[pheno * n_words];
      for (std::size_t r_id = n_loc_tot_ - n_loc_pheno; r_id < n_loc_tot_; r_id++) {
        std::size_t l_id = rng_unf_.sample(r_id + 1);
        std::size_t lw = l_id / 64; std::size_t lo = l_id % 64;
        std::size_t rw = r_id / 64; std::size_t ro = r_id % 64;
        if (loc_ptr[lw] & (1ull << lo)) loc_ptr[rw] |= (1ull << ro);
        else loc_ptr[lw] |= (1ull << lo);
      }
    }

    // optimise the expected panmictic genetic correlation matrix using greedy
    // approach
  }
}
