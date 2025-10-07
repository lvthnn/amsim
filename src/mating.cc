#include <cstddef>
#include <vector>
#include <stdexcept>
#include <numeric>
#include <algorithm>

#if defined(__APPLE__) && defined(USE_BLAS)
  #include <Accelerate/Accelerate.h>
#elif defined(__linux__) && defined(USE_BLAS)
  #include <cblas.h>
#endif

#include <Accelerate/Accelerate.h>

#include <amsim/mating.h>

namespace amsim::mating {
  std::vector<std::size_t> MatingModel::rand_state_() {
    std::vector<std::size_t> state_(n_sex_);
    std::iota(state_.begin(), state_.end(), n_sex_ + 1);
    std::shuffle(state_.begin(), state_.end(), g_);
    return state_;
  }

  GeneralModel::GeneralModel(const std::vector<double*> &ptr_tot,
                             std::vector<double> cor, const std::size_t n_itr,
                             const std::size_t n_sex, double tmp_init,
                             double tmp_decay)
    : MatingModel(MatingType::ASSORTATIVE, n_sex),
      ptr_tot_(ptr_tot),
      cor_(std::move(cor)),
      n_pheno_(ptr_tot.size()),
      n_itr_(n_itr),
      n_sex_(n_sex),
      tmp_init_(tmp_init),
      tmp_decay_(tmp_decay) {
    male_.resize(n_sex_ * n_pheno_);
    female_.resize(n_sex_ * n_pheno_);
  }

  // pack male and (permuted) female phenotypes into allocated buffers
  // male_ and female_
  void setup_(std::vector<std::size_t> state) {
    for (std::size_t pheno = 0; pheno < n_pheno_; pheno++) {
      const double* ptr_male_ = ptr_tot_[pheno];
      const double* ptr_female_ = ptr_male_ + n_sex_;
      double* male_col_ = male_ + pheno * n_sex_;
      double* female_col_ = female_ + pheno * n_sex_;

      std::copy(ptr_male_, ptr_male_ + n_sex_, male_col_);

      for (std::size_t ind = 0; ind < n_sex; ind++)
        female_col_[ind] = ptr_female_[state[ind]];
    }
  }

  std::vector<double> GeneralModel::cmp_cor_(std::vector<std::size_t> state) {
    std::vector<double> cross_cor_(n_pheno_ * n_pheno_);
    #if defined(USE_BLAS)
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  n_sex_, n_sex_, n_pheno_, 1.0 / static_cast<double>(n_pheno_),
                  male_.data(), 1, female_.data(), 1, 0.0, cross_cor_.data())
    #else
      throw std::runtime_error("Not implemented.");
    #endif
  }
}
