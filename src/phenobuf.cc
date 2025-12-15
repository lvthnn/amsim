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
#include <amsim/stats.h>

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

  std::vector<double> mean_male(n_pheno_);
  std::vector<double> mean_female(n_pheno_);
  std::vector<double> std_male(n_pheno_);
  std::vector<double> std_female(n_pheno_);

  for (std::size_t pheno = 0; pheno < n_pheno_; pheno++) {
    const double* tot_male = (*this)(ComponentType::TOTAL) + (pheno * n_ind_);
    const double* tot_female = tot_male + n_sex;

    mean_male[pheno] = stats::mean(n_sex, tot_male, 1);
    std_male[pheno] = stats::std(n_sex, tot_male, 1);
    mean_female[pheno] = stats::mean(n_sex, tot_female, 1);
    std_female[pheno] = stats::std(n_sex, tot_female, 1);
  }

  // compute latent score for each component
  for (ComponentType type :
       {ComponentType::GENETIC,
        ComponentType::ENVIRONMENTAL,
        ComponentType::VERTICAL,
        ComponentType::TOTAL}) {
    const double* buf_male = (*this)(type);
    const double* buf_female = buf_male + n_sex;
    double* buf_lat = (*this).latent(type);

    std::vector<double> pheno_stdbuf_male(n_sex * n_pheno_);
    std::vector<double> pheno_stdbuf_female(n_sex * n_pheno_);

    for (std::size_t pheno = 0; pheno < n_pheno_; ++pheno) {
      const double* pheno_male = buf_male + (pheno * n_ind_);
      const double* pheno_female = buf_female + (pheno * n_ind_);
      double* pheno_std_male = pheno_stdbuf_male.data() + (pheno * n_sex);
      double* pheno_std_female = pheno_stdbuf_female.data() + (pheno * n_sex);

      stats::standardise(
          n_sex,
          pheno_male,
          1,
          pheno_std_male,
          1,
          mean_male[pheno],
          std_male[pheno]);

      stats::standardise(
          n_sex,
          pheno_female,
          1,
          pheno_std_female,
          1,
          mean_female[pheno],
          std_female[pheno]);
    }

    cblas_dgemm(
        CblasColMajor,
        CblasNoTrans,
        CblasNoTrans,
        n_sex,
        n_pheno_,
        n_pheno_,
        1.0,
        pheno_stdbuf_male.data(),
        n_sex,
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
        pheno_stdbuf_female.data(),
        n_sex,
        VT.data(),
        n_pheno_,
        0.0,
        buf_lat + n_sex,
        n_ind_);
  }
}

}  // namespace amsim
