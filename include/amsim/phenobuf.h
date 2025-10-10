#ifndef AMSIMCPP_PHENOBUF_H
#define AMSIMCPP_PHENOBUF_H

#include <cstddef>
#include <cstdint>

#include <amsim/componenttype.h>

namespace amsim {
  class PhenoBuf {
  public:
    PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno);

    inline const double* operator()(std::size_t id, ComponentType type) const {
      if (occupied_[id])
        throw std::runtime_error("buffer slot already occupied");
      switch (type) {
        case ComponentType::GENETIC: return &buffer_[id * n_ind_];
        case ComponentType::ENVIRONMENTAL: return &buffer_[(n_pheno_ + id) * n_ind_];
        case ComponentType::VERTICAL: return &buffer_[(2 * n_pheno_ + id) * n_ind_];
        case ComponentType::TOTAL: return &buffer_[(3 * n_pheno_ + id) * n_ind_];
      }
    }

    inline double* operator()(std::size_t id, ComponentType type) {
      if (occupied_[id])
        throw std::runtime_error("buffer slot already occupied");
      switch (type) {
        case ComponentType::GENETIC: return &buffer_[id * n_ind_];
        case ComponentType::ENVIRONMENTAL: return &buffer_[(n_pheno_ + id) * n_ind_];
        case ComponentType::VERTICAL: return &buffer_[(2 * n_pheno_ + id) * n_ind_];
        case ComponentType::TOTAL: return &buffer_[(3 * n_pheno_ + id) * n_ind_];
      }
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