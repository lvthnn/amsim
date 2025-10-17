#include <cstddef>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>

#if defined(__APPLE__) && defined(USE_BLAS)
  #include <Accelerate/Accelerate.h>
#elif defined(__linux__) && defined(USE_BLAS)
  #include <cblas.h>
#endif

#include <amsim/rng.h>
#include <amsim/mating.h>

namespace amsim {
  std::vector<std::size_t> MatingModel::rand_state_() {
    std::vector<std::size_t> state_(n_sex_);
    std::iota(state_.begin(), state_.end(), 0);
    std::shuffle(state_.begin(), state_.end(), g_);
    return state_;
  }

  AssortativeModel::AssortativeModel(const PhenotypeList &phenotypes,
                                     std::vector<double> cor,
                                     const std::size_t n_itr,
                                     const std::size_t n_sex,
                                     const rng::Xoshiro256ss &rng,
                                     double tmp_init, double tmp_decay)
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
      state_(rand_state_()),
      swap_(rng),
      acc_(rng) {
    male_.resize(n_sex_ * n_pheno_);
    female_.resize(n_sex_ * n_pheno_);
  }

  void AssortativeModel::arrange_() {
    std::vector<double> ones(n_sex_, 1.0);
    for (std::size_t pheno = 0; pheno < n_pheno_; pheno++) {
      const double* ptr_m = ptr_tot_[pheno];
      const double* ptr_f = ptr_m + n_sex_;

      const double scale   = 1.0 / static_cast<double>(n_sex_);
      const double mean_m  = scale * cblas_ddot(n_sex_, ptr_m, 1, ones.data(), 1);
      const double mean_f  = scale * cblas_ddot(n_sex_, ptr_f, 1, ones.data(), 1);
      const double sumsq_m = scale * cblas_ddot(n_sex_, ptr_m, 1, ptr_m, 1);
      const double sumsq_f = scale * cblas_ddot(n_sex_, ptr_f, 1, ptr_f, 1);
      const double sd_m    = std::sqrt(sumsq_m - mean_m * mean_m);
      const double sd_f    = std::sqrt(sumsq_f - mean_f * mean_f);
      // compute the sex-segregated means and standard deviations
      // compute the sex-segregated means and standard deviations

      for (std::size_t ind = 0; ind < n_sex_; ind++) {
        male_[pheno * n_sex_ + ind]   = (ptr_m[ind] - mean_m) / sd_m;
        female_[pheno * n_sex_ + ind] = (ptr_f[state_[ind]] - mean_f) / sd_f; 
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
      for (std::size_t p2 = 0; p2 < n_pheno_; p2++) {
        std::size_t pair = p1 * n_pheno_ + p2;
        double m0 = male_[p1 * n_sex_ + i0];
        double m1 = male_[p1 * n_sex_ + i1];
        double f0 = female_[p2 * n_sex_ + i0];
        double f1 = female_[p2 * n_sex_ + i1];
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
    for (std::size_t el = 0; el < cor_mat.size(); el++) {
      if (el % n_pheno_ == 0) std::cout << "\n";
      std::cout << cor_mat[el] << "\t";
    }
    std::cout << "\n";
  }

  std::vector<std::size_t> AssortativeModel::match() {
    if (n_itr_ == 0) return state_;
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
        std::swap(state_[i0], state_[i1]);

        for (std::size_t pheno = 0; pheno < n_pheno_; pheno++)
          std::swap(female_[pheno * n_sex_ + i0], female_[pheno * n_sex_ + i1]);

        cblas_daxpy(dim, 1.0, delta.data(), 1, cur.data(), 1);
      }

      tmp_cur *= tmp_decay_;
    }

    return state_;
  }

  void AssortativeModel::update(const PhenotypeList &phenotypes) {
    std::vector<const double*> ptr_new;
    for (const auto &pheno : phenotypes)
      ptr_new.push_back(pheno(ComponentType::TOTAL));
    ptr_tot_ = ptr_new;
  }
}
