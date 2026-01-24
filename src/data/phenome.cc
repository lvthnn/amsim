#include <amsim/data/phenome.h>

#include <cstddef>

namespace amsim::phenome {

PhenoBuf::PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno)
    : n_ind_(n_ind), n_sex_(n_ind_ / 2), n_pheno_(n_pheno) {
  data_.resize(4 * n_ind_ * n_pheno_);
}

void PhenoBuf::compute_stats() {
  for (ComponentType type :
       {ComponentType::GENETIC,
        ComponentType::ENVIRONMENTAL,
        ComponentType::VERTICAL,
        ComponentType::TOTAL}) {
    comp_mean_ = (*this)(type).colwise().mean();
    comp_var_ =
        ((*this)(type).rowwise() - comp_mean_).array().square().colwise().sum();
    comp_var_ *= (1.0 / n_ind_);
  }
}

}  // namespace amsim::phenome
