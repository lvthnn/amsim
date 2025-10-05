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

  GeneralModel::GeneralModel(std::vector<std::vector<double> const*> vals_ptr,
                             std::vector<double> cor, const std::size_t n_itr,
                             const std::size_t n_sex, double tmp_init,
                             double tmp_decay)
    : MatingModel(MatingType::ASSORTATIVE, n_sex),
      cor_(std::move(cor)),
      vals_ptr_(std::move(vals_ptr)),
      n_itr_(n_itr),
      tmp_init_(tmp_init),
      tmp_decay_(tmp_decay) {
    if (cor_.size() != std::pow(vals_ptr_.size(), 2))
      throw std::runtime_error("cor_ length must be vals_ptr_.size() squared");
  }

  std::vector<double> GeneralModel::cmp_cor_(std::vector<std::size_t> state) {
    std::size_t dim = vals_ptr_.size();
    std::size_t n_el = dim * dim;
    std::vector<double> cor(n_el, 0.0);

    for (std::size_t el = 0; el < n_el; el++) {
      std::size_t r = el / dim;
      std::size_t c = el % dim;

      #if defined(USE_BLAS)
        cor[el] = cblas_ddot(2 * n_sex_, (*vals_ptr_[r]).data(), 1, (*vals_ptr_[c]).data(), 1);
      #else
        cor[el] = ddot_(2 * n_sex_, double *dx, int *incx, double *dy, int *incy)
      #endif
    }

    return cor;
  }
}
