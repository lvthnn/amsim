#ifndef AMSIMCPP_PHENOBUF_H
#define AMSIMCPP_PHENOBUF_H

#include <cstddef>
#include <vector>

#include <amsim/componenttype.h>

namespace amsim {
  class PhenoBuf {
  public:
    PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno);

    inline const double* operator()(std::size_t id, ComponentType type) const {
      return &buffer_[n_ind_ * static_cast<int>(type) * n_pheno_ + id];
    }

    inline double* operator()(std::size_t id, ComponentType type) {
      return &buffer_[n_ind_ * static_cast<int>(type) * n_pheno_ + id];
    }

    inline std::size_t n_ind() const noexcept { return n_ind_; }
    inline bool occupied(std::size_t id) const noexcept { return occupied_[id]; }

    std::optional<std::size_t> unoccupied() const;
    void occupy(std::size_t);

  private:
    const std::size_t n_ind_;
    const std::size_t n_pheno_;
    std::vector<double> buffer_;
    std::vector<bool> occupied_;
  };
}

#endif // AMSIMCPP_PHENOBUF_H
