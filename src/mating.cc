#include <cstddef>
#include <vector>
#include <stdexcept>
#include <numeric>
#include <algorithm>

#include <amsim/mating.h>

namespace amsim::mating {
  std::vector<std::size_t> MatingModel::rand_state_() {
    std::vector<std::size_t> state_(n_sex_);
    std::iota(state_.begin(), state_.end(), n_sex_ + 1);
    std::shuffle(state_.begin(), state_.end(), g_);
    return state_;
  }

  GeneralModel::GeneralModel(std::vector<std::vector<double> const*> vals_ptr,
                             std::tuple<std::size_t, std::size_t, double> cor,
                             const std::size_t n_itr, const std::size_t n_sex,
                             double tmp_init, double tmp_decay)
    : MatingModel(MatingType::ASSORTATIVE, n_sex),
      vals_ptr_(std::move(vals_ptr)),
      cor_(cor),
      n_itr_(n_itr),
      n_sex_(n_sex),
      tmp_init_(tmp_init),
      tmp_decay_(tmp_decay) {}

  double GeneralModel::delta_(std::vector<std::size_t> cur, std::size_t i0,
                              std::size_t i1) {
    // compute the delta vector
    double delta = 0.0;
    for (std::size_t el = 0; el < vals_ptr_.size(); el++) {
      double m0 = (*vals_ptr_[el])[i0]; double m1 = (*vals_ptr_[el])[i1];
      double df = (*vals_ptr_[el])[cur[i1]] - (*vals_ptr_[el])[cur[i0]];
      delta += m0 * df - m1 * df;
    }

    // get the energy differential using dot products

    return delta;
  }

  void GeneralModel::update_vals(std::vector<std::vector<double> const*> vals_ptr) {
    if (vals_ptr.size() != vals_ptr_.size())
      throw std::runtime_error("vals_ptr must be same length as vals_ptr_");

    for (std::size_t el = 0; el < vals_ptr_.size(); el++)
      vals_ptr_[el] = vals_ptr[el];
  }

  std::vector<std::size_t> GeneralModel::match() {
    std::vector<std::size_t> state = rand_state_();
    std::uniform_int_distribution<std::size_t> unif(0, n_sex_ - 1); 
    double tmp = tmp_init_;

    for (std::size_t it = 0; it < n_itr_; it++) {
      tmp *= tmp_decay_;
    }
    return state;
  }
}
