#ifndef AMSIMCPP_PHENOBUF_H
#define AMSIMCPP_PHENOBUF_H

#include <amsim/component_type.h>

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

namespace amsim {
class PhenoBuf {
 public:
  PhenoBuf(std::size_t n_ind, std::size_t n_pheno, bool require_lat);

  const double* operator()(std::size_t id, ComponentType type) const {
    return &buf_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
  }

  double* operator()(std::size_t id, ComponentType type) {
    return &buf_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
  }

  double* operator()(ComponentType type) {
    return &buf_[n_ind_ * n_pheno_ * static_cast<int>(type)];
  }

  const double* latent(std::size_t id, ComponentType type) const {
    if (buf_lat_.empty()) throw std::runtime_error("latent buffer not in use");
    return &buf_lat_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
  }

  double* latent(ComponentType type) {
    if (buf_lat_.empty()) throw std::runtime_error("latent buffer not in use");
    return &buf_lat_[n_ind_ * n_pheno_ * static_cast<int>(type)];
  }

  std::size_t n_ind() const noexcept { return n_ind_; }
  bool occupied(std::size_t id) const noexcept { return occupied_[id]; }
  bool has_lat() const noexcept { return !buf_lat_.empty(); }

  std::optional<std::size_t> unoccupied() const;
  void occupy(std::size_t);
  void score_latent(
      const std::vector<double>& U, const std::vector<double>& VT);

 private:
  const std::size_t n_ind_;
  const std::size_t n_pheno_;
  std::vector<double> buf_;
  std::vector<double> buf_lat_;
  std::vector<bool> occupied_;
};

}  // namespace amsim

#endif  // AMSIMCPP_PHENOBUF_H
