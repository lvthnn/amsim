#include <algorithm>
#include <cstddef>
#include <optional>
#include <vector>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#include <amsim/phenobuf.h>

namespace amsim {
PhenoBuf::PhenoBuf(
    const std::size_t n_ind, const std::size_t n_pheno, bool require_lat)
    : n_ind_(n_ind), n_pheno_(n_pheno) {
  buf_.resize(4 * n_ind_ * n_pheno_);
  if (require_lat) buf_lat_.resize(4 * n_ind_ * n_pheno_);
  occupied_.resize(n_pheno_, false);
}

std::optional<std::size_t> PhenoBuf::unoccupied() const {
  auto it = std::ranges::find(occupied_, false);
  if (it != occupied_.end()) {
    std::size_t id = std::distance(occupied_.begin(), it);
    return id;
  }
  return std::nullopt;
}

void PhenoBuf::occupy(std::size_t id) { occupied_[id] = true; }

void PhenoBuf::score_latent(
    const std::vector<double>& U, const std::vector<double>& VT) {
  const std::size_t n_sex = n_ind_ / 2;
  const int lda_buf = 2 * n_sex;

  // compute latent score for each component
  for (ComponentType type :
       {ComponentType::GENETIC,
        ComponentType::ENVIRONMENTAL,
        ComponentType::VERTICAL,
        ComponentType::TOTAL}) {
    const double* buf_male = (*this)(type);
    const double* buf_female = buf_male + n_sex;
    double* buf_lat = (*this).latent(type);

    cblas_dgemm(
        CblasColMajor,
        CblasNoTrans,
        CblasNoTrans,
        n_sex,
        n_pheno_,
        n_pheno_,
        1.0,
        buf_male,
        lda_buf,
        U.data(),
        n_pheno_,
        0.0,
        buf_lat,
        n_ind_);

    cblas_dgemm(
        CblasColMajor,
        CblasNoTrans,
        CblasTrans,
        n_sex,
        n_pheno_,
        n_pheno_,
        1.0,
        buf_female,
        lda_buf,
        VT.data(),
        n_pheno_,
        0.0,
        buf_lat + n_sex,
        n_ind_);
  }
}
}  // namespace amsim
