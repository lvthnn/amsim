#ifndef AMSIMCPP_PHENOBUF_H
#define AMSIMCPP_PHENOBUF_H

#pragma once
#include <cstddef>
#include <vector>

#include <amsim/component_type.h>

namespace amsim {
  class PhenoBuf {
  public:
    PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno, bool require_lat);

    inline const double* operator()(std::size_t id, ComponentType type) const {
      return &buf_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    }

    inline double* operator()(std::size_t id, ComponentType type) {
      return &buf_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    }

    inline double* operator()(ComponentType type) {
      return &buf_[n_ind_ * n_pheno_ * static_cast<int>(type)];
    }

    inline const double* latent(std::size_t id, ComponentType type) const {
      if (buf_lat_.empty()) throw std::runtime_error("latent buffer not in use");
      return &buf_lat_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    }

    inline std::size_t n_ind() const noexcept { return n_ind_; }
    inline bool occupied(std::size_t id) const noexcept { return occupied_[id]; }

    std::optional<std::size_t> unoccupied() const;
    void occupy(std::size_t);
    void score_latent(const std::vector<double>& U, const std::vector<double>& VT);
    const bool require_lat;

  private:
    const std::size_t n_ind_;
    const std::size_t n_pheno_;
    std::vector<double> buf_;
    std::vector<double> buf_lat_;
    std::vector<bool> occupied_;
  };
}

#endif // AMSIMCPP_PHENOBUF_H
