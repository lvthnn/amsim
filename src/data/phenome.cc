#include <amsim/data/phenome.h>

#include <cstddef>

namespace amsim::phenome {

PhenoBuf::PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno)
    : n_ind_(n_ind),
      n_sex_(n_ind_ / 2),
      n_pheno_(n_pheno),
      data_(4 * n_ind_ * n_pheno_),
      comp_mean_(4, n_pheno_),
      comp_var_(4, n_pheno_) {};

void PhenoBuf::compute_stats() {
  for (ComponentType type :
       {ComponentType::GENETIC,
        ComponentType::ENVIRONMENTAL,
        ComponentType::VERTICAL,
        ComponentType::TOTAL}) {
    Eigen::Index col_type = static_cast<int>(type);
    comp_mean_.row(col_type) = (*this)(type).colwise().mean();
    comp_var_.row(col_type) =
        ((*this)(type).rowwise() - comp_mean_.row(col_type))
            .array()
            .square()
            .colwise()
            .sum() /
        static_cast<double>(n_ind_);
  }
}

}  // namespace amsim::phenome
