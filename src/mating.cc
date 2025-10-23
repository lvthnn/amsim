#include <cstddef>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <iomanip>

#ifdef __APPLE__
  #include <Accelerate/Accelerate.h>
#else
  #include <cblas.h>
  #include <lapacke.h>
#endif

#include <amsim/rng.h>
#include <amsim/utils.h>
#include <amsim/mating.h>

namespace amsim {
  std::vector<std::size_t> MatingModel::rand_state_() {
    std::vector<std::size_t> state(n_sex_);
    std::iota(state.begin(), state.end(), 0);
    std::shuffle(state.begin(), state.end(), g_);
    return state;
  }

  AssortativeModel::AssortativeModel(const PhenotypeList&     phenotypes,
                                     std::vector<double>      cor,
                                     const std::size_t        n_itr,
                                     const std::size_t        n_sex,
                                     const rng::Xoshiro256ss& rng,
                                     double                   tmp_init,
                                     double                   tmp_decay)
    : MatingModel(MatingType::ASSORTATIVE, n_sex),
      ptr_tot_([&]() {
        std::vector<const double*> res;
        res.reserve(phenotypes.size());
        for (const auto &pheno : phenotypes)
          res.push_back(pheno(ComponentType::TOTAL));
        return res;
      }()),
      cor_(std::move(cor)),
      n_pheno_(phenotypes.size()),
      n_sex_(n_sex),
      n_itr_(n_itr),
      tmp_init_(tmp_init),
      tmp_decay_(tmp_decay),
      fuzz_(rng),
      swap_(rng),
      acc_(rng),
      cor_S(n_pheno_),
      cor_U(n_pheno_ * n_pheno_),
      cor_VT(n_pheno_ * n_pheno_),
      state(n_sex_) {
    // copy since dgesvd from LAPACK destroys original matrix
    std::vector<double> cor_copy_ = cor_;

    char   clpk_job      = 'A';
    int    clpk_n_pheno_ = static_cast<int>(n_pheno_);
    int    clpk_lwork    = -1;
    int    clpk_info;
    double clpk_wkopt;

    // perform workspace query
    #if defined(__APPLE__)
      dgesvd_(&clpk_job, &clpk_job, &clpk_n_pheno_, &clpk_n_pheno_, cor_copy_.data(),
              &clpk_n_pheno_, cor_S.data(), cor_U.data(), &clpk_n_pheno_,
              cor_VT.data(), &clpk_n_pheno_, &clpk_wkopt, &clpk_lwork, &clpk_info);
    #else
      LAPACK_dgesvd(&clpk_job, &clpk_job, &clpk_n_pheno_, &clpk_n_pheno_, cor_copy_.data(),
                    &clpk_n_pheno_, cor_S.data(), cor_U.data(), &clpk_n_pheno_,
                    cor_VT.data(), &clpk_n_pheno_, &clpk_wkopt, &clpk_lwork, &clpk_info);
    #endif

    // allocate workspace buffer
    clpk_lwork = static_cast<int>(clpk_wkopt);
    std::vector<double> clpk_work(clpk_lwork);

    // perform SDV computation
    #if defined(__APPLE__)
      dgesvd_(&clpk_job, &clpk_job, &clpk_n_pheno_, &clpk_n_pheno_, cor_copy_.data(),
              &clpk_n_pheno_, cor_S.data(), cor_U.data(), &clpk_n_pheno_,
              cor_VT.data(), &clpk_n_pheno_, clpk_work.data(), &clpk_lwork, &clpk_info);
    #else
      LAPACK_dgesvd(&clpk_job, &clpk_job, &clpk_n_pheno_, &clpk_n_pheno_, cor_copy_.data(),
                    &clpk_n_pheno_, cor_S.data(), cor_U.data(), &clpk_n_pheno_,
                    cor_VT.data(), &clpk_n_pheno_, clpk_work.data(), &clpk_lwork, &clpk_info);
    #endif

    male_.resize(n_sex_ * n_pheno_);
    female_.resize(n_sex_ * n_pheno_);
  }

  void AssortativeModel::arrange_() {
    std::vector<double> ones(n_sex_, 1.0);
    for (std::size_t pheno = 0; pheno < n_pheno_; pheno++) {
      const double* ptr_m = ptr_tot_[pheno];
      const double* ptr_f = ptr_m + n_sex_;

      // @TODO: replace this with calls to stats header
      const double scale   = 1.0 / static_cast<double>(n_sex_);
      const double mean_m  = scale * cblas_ddot(n_sex_, ptr_m, 1, ones.data(), 1);
      const double mean_f  = scale * cblas_ddot(n_sex_, ptr_f, 1, ones.data(), 1);
      const double sumsq_m = scale * cblas_ddot(n_sex_, ptr_m, 1, ptr_m, 1);
      const double sumsq_f = scale * cblas_ddot(n_sex_, ptr_f, 1, ptr_f, 1);
      const double sd_m    = std::sqrt(sumsq_m - mean_m * mean_m);
      const double sd_f    = std::sqrt(sumsq_f - mean_f * mean_f);

      for (std::size_t ind = 0; ind < n_sex_; ind++) {
        male_[pheno * n_sex_ + ind]   = (ptr_m[ind] - mean_m) / sd_m;
        female_[pheno * n_sex_ + ind] = (ptr_f[state[ind]] - mean_f) / sd_f;
      }
    }
  }

  std::vector<double> AssortativeModel::compute_cor_() {
    std::vector<double> cross_cor_(n_pheno_ * n_pheno_);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                n_pheno_, n_pheno_, n_sex_, 1.0 / static_cast<double>(n_sex_),
                male_.data(), n_sex_, female_.data(), n_sex_, 0.0, cross_cor_.data(), n_pheno_);
    return cross_cor_;
  }

  std::vector<double> AssortativeModel::compute_delta_(std::size_t i0, std::size_t i1) {
    std::vector<double> res(n_pheno_ * n_pheno_);
    const double scale = 1.0 / static_cast<double>(n_sex_);

    for (std::size_t p1 = 0; p1 < n_pheno_; p1++) {
      double m0 = male_[p1 * n_sex_ + i0];
      double m1 = male_[p1 * n_sex_ + i1];
      for (std::size_t p2 = 0; p2 < n_pheno_; p2++) {
        double f0 = female_[p2 * n_sex_ + i0];
        double f1 = female_[p2 * n_sex_ + i1];
        // column-major: index = row + col * n_rows
        std::size_t pair = p1 + p2 * n_pheno_;
        res[pair] = scale * (m0 * (f1 - f0) + m1 * (f0 - f1));
      }
    }
    return res;
  }

  double AssortativeModel::compute_denergy_(const std::vector<double> &cur,
                                            const std::vector<double> &target,
                                            const std::vector<double> &delta) {
    std::size_t dim = n_pheno_ * n_pheno_;
    std::vector<double> diff = cur;

    cblas_daxpy(dim, -1.0, target.data(), 1, diff.data(), 1);
    return cblas_ddot(dim, delta.data(), 1, delta.data(), 1) +
           2.0 * cblas_ddot(dim, diff.data(), 1, delta.data(), 1);
  }

  void AssortativeModel::display_cor() {
    std::vector<double> cor_mat = compute_cor_();
    std::cerr << std::setprecision(5);
    for (std::size_t el = 0; el < cor_mat.size(); el++) {
      if (el % n_pheno_ == 0) std::cerr << "\n";
      std::cerr << cor_mat[el] << "\t";
    }
    std::cerr << "\n";
  }

  void AssortativeModel::init_state() {
    // compute dominant latent phenotypes
    double latent_cor   = cor_S[0];
    double latent_noise = std::sqrt(1.0 - latent_cor * latent_cor);
    std::vector<double> latent_male(n_sex_);
    std::vector<double> latent_female(n_sex_);

    int lda_tot = 2 * n_sex_;

    cblas_dgemv(CblasColMajor, CblasNoTrans, n_sex_, n_pheno_, 1.0, ptr_tot_[0], lda_tot,
                cor_U.data(), 1, 0.0, latent_male.data(), 1);

    cblas_dgemv(CblasColMajor, CblasNoTrans, n_sex_, n_pheno_, 1.0, ptr_tot_[0] + n_sex_,
                lda_tot, cor_VT.data(), n_pheno_, 0.0, latent_female.data(), 1);

    // add noise to get correlation to equal latent_cor
    std::vector<double> latent_fuzz(n_sex_);
    fuzz_.fill(latent_fuzz.data(), n_sex_);
    for (std::size_t ind = 0; ind < n_sex_; ind++)
      latent_female[ind] += latent_noise * latent_fuzz[ind];

    // generate initial state based on matching of order statistics
    std::vector<std::size_t> idx_male   = utils::order(latent_male);
    std::vector<std::size_t> idx_female = utils::order(latent_female);

    // assemble the initial state
    for (std::size_t rank = 0; rank < n_sex_; ++rank)
      state[idx_male[rank]] = idx_female[rank];
  }

  std::vector<std::size_t> AssortativeModel::match() {
    if (n_itr_ == 0) return state;
    const std::size_t dim = n_pheno_ * n_pheno_;
    double tmp_cur = tmp_init_;

    arrange_();
    std::vector<double> cur = compute_cor_();

    for (std::size_t itr = 0; itr < n_itr_; itr++) {
      std::size_t i0 = swap_.sample(n_sex_);
      std::size_t i1 = swap_.sample(n_sex_);
      while (i0 == i1) i1 = swap_.sample(n_sex_);

      std::vector<double> delta = compute_delta_(i0, i1);

      double denergy  = compute_denergy_(cur, cor_, delta);
      double acc_prob = std::min(1.0, std::exp(-denergy / tmp_cur));
      double u        = acc_.sample(1.0);

      if (u < acc_prob) {
        std::swap(state[i0], state[i1]);

        for (std::size_t pheno = 0; pheno < n_pheno_; pheno++)
          std::swap(female_[pheno * n_sex_ + i0], female_[pheno * n_sex_ + i1]);

        cblas_daxpy(dim, 1.0, delta.data(), 1, cur.data(), 1);
      }

      tmp_cur *= tmp_decay_;
    }

    return state;
  }

  void AssortativeModel::update(const PhenotypeList &phenotypes) {
    std::vector<const double*> ptr_new;
    for (const auto &pheno : phenotypes)
      ptr_new.push_back(pheno(ComponentType::TOTAL));
    ptr_tot_ = ptr_new;
  }
}
