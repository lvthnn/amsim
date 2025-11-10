#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

#include <amsim/phenoarch.h>
#include <amsim/rng.h>

namespace amsim {
PhenoArch::PhenoArch(
    std::size_t n_pheno,
    std::size_t n_loc_total,
    std::vector<std::size_t> n_loc,
    std::vector<double> h2_gen,
    std::vector<double> h2_env,
    std::vector<double> gen_cor,
    std::vector<double> env_cor,
    const rng::Xoshiro256ss& rng)
    : n_pheno_(n_pheno),
      n_loc_tot_(n_loc_total),
      n_loc_(std::move(n_loc)),
      h2_gen_(std::move(h2_gen)),
      h2_env_(std::move(h2_env)),
      gen_cor_(std::move(gen_cor)),
      env_chol_(std::move(env_cor)),
      rng_polar_(rng),
      rng_unf_(rng) {
  assert(n_loc_.size() == n_pheno_);
  assert(h2_gen_.size() == n_pheno_);
  assert(gen_cor_.size() == n_pheno_ * n_pheno_);
  assert(env_chol_.size() == n_pheno_ * n_pheno_);

  // cast to LAPACK-legible form
  char clpk_uplo_ = 'L';
  int clpk_n_pheno_ = static_cast<int>(n_pheno_);
  int clpk_lda_ = static_cast<int>(n_pheno_);
  int clpk_out_;

#if defined(__APPLE__)
  dpotrf_(
      &clpk_uplo_, &clpk_n_pheno_, env_chol_.data(), &clpk_lda_, &clpk_out_);
#else
  LAPACKE_dpotrf(
      &clpk_uplo_, &clpk_n_pheno_, env_chol_.data(), &clpk_lda_, &clpk_out_);
#endif

  for (std::size_t c = 0; c < n_pheno; ++c)
    for (std::size_t r = 0; r < c; ++r) env_chol_[c * n_pheno_ + r] = 0.0;

  // add initial state optimisation to make this more reliable
  optim_arch(1e4);
}

void PhenoArch::gen_env(double* ptr_env, const std::size_t n_ind) {
  rng_polar_.fill(ptr_env, n_ind * n_pheno_);
  cblas_dtrmm(
      CblasColMajor,
      CblasRight,
      CblasLower,
      CblasTrans,
      CblasNonUnit,
      n_ind,
      n_pheno_,
      1.0,
      env_chol_.data(),
      n_pheno_,
      ptr_env,
      n_ind);

  // scale environmental components along margins
  for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno)
    cblas_dscal(n_ind, std::sqrt(h2_env_[pheno]), &ptr_env[pheno * n_ind], 1);
}

void PhenoArch::print_correlations(
    const std::vector<std::size_t>& intersect) const {
  std::cerr << "Phenotype correlations:\n";
  std::cerr << std::fixed << std::setprecision(4);

  std::size_t id = 0;
  for (std::size_t i = 0; i < n_pheno_; ++i) {
    for (std::size_t j = i + 1; j < n_pheno_; ++j) {
      const double actual_cor =
          intersect[id] / std::sqrt(n_loc_[i] * n_loc_[j]);
      std::cerr << "  Pheno " << i << " vs " << j << ": " << actual_cor << "\n";
      ++id;
    }
  }
  std::cerr << "\n";
}

std::vector<std::uint64_t> PhenoArch::init_mask_() {
  const std::size_t n_words = (n_loc_tot_ + 63) / 64;
  std::vector<std::uint64_t> loc_mask(n_words * n_pheno_);

  for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno) {
    std::uint64_t* loc_ptr = &loc_mask[pheno * n_words];

    const std::size_t n_loc_pheno = n_loc_[pheno];
    for (std::size_t r_id = n_loc_tot_ - n_loc_pheno; r_id < n_loc_tot_;
         ++r_id) {
      const std::size_t l_id = rng_unf_.sample(r_id + 1);
      const std::size_t lw = l_id / 64;
      const std::size_t lo = l_id % 64;
      const std::size_t rw = r_id / 64;
      const std::size_t ro = r_id % 64;

      if (loc_ptr[lw] & (1ull << lo))
        loc_ptr[rw] |= (1ull << ro);
      else
        loc_ptr[lw] |= (1ull << lo);
    }
  }
  return loc_mask;
}

std::vector<std::size_t> PhenoArch::init_intersect_(
    const std::vector<std::uint64_t>& mask) const {
  std::vector<std::size_t> intersect(n_pheno_ * (n_pheno_ - 1) / 2);
  const std::size_t n_words = (n_loc_tot_ + 63) / 64;

  std::size_t id = 0;
  for (std::size_t i = 0; i < n_pheno_; ++i) {
    for (std::size_t j = i + 1; j < n_pheno_; ++j) {
      const std::uint64_t* pheno_i = &mask[n_words * i];
      const std::uint64_t* pheno_j = &mask[n_words * j];
      for (std::size_t word = 0; word < n_words; ++word) {
        intersect[id] += __builtin_popcountll(pheno_i[word] & pheno_j[word]);
      }
      ++id;
    }
  }
  return intersect;
}

std::vector<double> PhenoArch::init_weights_() const {
  std::vector<double> weights((n_pheno_ * (n_pheno_ - 1)) / 2);
  std::size_t id = 0;
  for (std::size_t i = 0; i < n_pheno_; ++i) {
    for (std::size_t j = i + 1; j < n_pheno_; ++j) {
      weights[id] =
          gen_cor_[j * n_pheno_ + i] * std::sqrt(n_loc_[i] * n_loc_[j]);
      ++id;
    }
  }
  return weights;
}

void PhenoArch::optim_arch(std::size_t max_it, double eps) {
  // @TODO: Make optimisation routine terminate for error < eps
  loc_mask_ = init_mask_();
  std::vector<std::size_t> intersect = init_intersect_(loc_mask_);
  std::vector<double> weights = init_weights_();

  const std::size_t n_words = (n_loc_tot_ + 63) / 64;
  for (std::size_t it = 0; it < max_it; ++it) {
    std::size_t pheno = it % n_pheno_;
    std::size_t opt_add = 0;
    std::size_t opt_del = 0;
    double opt_add_delta = std::numeric_limits<double>::max();
    double opt_del_delta = std::numeric_limits<double>::max();

    std::uint64_t* ptr_pheno = &loc_mask_[n_words * pheno];

    // scan through loci and determine optimal addition and deletion
    for (std::size_t loc = 0; loc < n_loc_tot_; ++loc) {
      std::size_t block = loc / 64;
      std::size_t offset = loc % 64;

      double delta = 0.0;
      bool causal = ptr_pheno[block] & (1ull << offset);

      for (std::size_t pheno_adj = 0; pheno_adj < n_pheno_; ++pheno_adj) {
        if (pheno == pheno_adj) continue;

        std::uint64_t* ptr_pheno_adj = &loc_mask_[n_words * pheno_adj];
        bool causal_adj = (ptr_pheno_adj[block] & (1ull << offset));

        if (causal_adj) {
          // convert triangular indices to linear index
          std::size_t i = std::min(pheno, pheno_adj);
          std::size_t j = std::max(pheno, pheno_adj);
          std::size_t idx = i * n_pheno_ - (i * (i + 1)) / 2 + (j - i - 1);
          std::size_t intersect_prev = intersect[idx];
          std::size_t intersect_cur = intersect_prev + (causal ? -1 : 1);
          double target = weights[idx];

          delta += (intersect_cur - target) * (intersect_cur - target) -
                   (intersect_prev - target) * (intersect_prev - target);
        }
      }

      if (causal && delta < opt_del_delta) {
        opt_del = loc;
        opt_del_delta = delta;
      } else if (!causal && delta < opt_add_delta) {
        opt_add = loc;
        opt_add_delta = delta;
      }
    }

    // apply the swap if it improves the cost
    double total_delta = opt_add_delta + opt_del_delta;
    if (total_delta < -eps) {
      // perform deletion
      std::size_t del_block = opt_del / 64;
      std::size_t del_offset = opt_del % 64;
      ptr_pheno[del_block] &= ~(1ull << del_offset);

      // perform addition
      std::size_t add_block = opt_add / 64;
      std::size_t add_offset = opt_add % 64;
      ptr_pheno[add_block] |= (1ull << add_offset);

      // update intersections and weights for deletion
      for (std::size_t pheno_adj = 0; pheno_adj < n_pheno_; ++pheno_adj) {
        if (pheno == pheno_adj) continue;

        std::uint64_t* ptr_pheno_adj = &loc_mask_[n_words * pheno_adj];
        if (ptr_pheno_adj[del_block] & (1ull << del_offset)) {
          std::size_t i = std::min(pheno, pheno_adj);
          std::size_t j = std::max(pheno, pheno_adj);
          std::size_t idx = i * n_pheno_ - (i * (i + 1)) / 2 + (j - i - 1);
          intersect[idx]--;
        }
      }

      // update intersections and weights for addition
      for (std::size_t pheno_adj = 0; pheno_adj < n_pheno_; ++pheno_adj) {
        if (pheno == pheno_adj) continue;

        std::uint64_t* ptr_pheno_adj = &loc_mask_[n_words * pheno_adj];
        if (ptr_pheno_adj[add_block] & (1ull << add_offset)) {
          std::size_t i = std::min(pheno, pheno_adj);
          std::size_t j = std::max(pheno, pheno_adj);
          std::size_t idx = i * n_pheno_ - (i * (i + 1)) / 2 + (j - i - 1);
          intersect[idx]++;
        }
      }
    }
  }
}

std::vector<std::size_t> PhenoArch::pheno_mask(
    const std::size_t pheno_id) const {
  const std::size_t n_words = (n_loc_tot_ + 63) / 64;
  const std::uint64_t* ptr_pheno = &loc_mask_[pheno_id * n_words];
  std::vector<std::size_t> loci;

  for (std::size_t word = 0; word < n_words; ++word) {
    for (std::uint64_t mask = ptr_pheno[word]; mask; mask &= (mask - 1)) {
      const std::size_t offset =
          static_cast<std::uint64_t>(__builtin_ctzll(mask));
      const std::size_t loc = 64 * word + offset;
      loci.push_back(loc);
    }
  }
  return loci;
}
}  // namespace amsim
