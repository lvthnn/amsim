//------------------------------------------------------------------------------
// amsimcpp : phenobuf.h
//------------------------------------------------------------------------------

#ifndef AMSIMCPP_PHENOBUF_H
#define AMSIMCPP_PHENOBUF_H

#include <cstddef>
#include <cstdint>

#include <amsim/componenttype.h>

class PhenoBuf {
public:
  PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno);

  inline const double* operator()(std::size_t id, ComponentType type) const noexcept {
    switch (type) {
      case ComponentType::GENETIC:
        return ptr_gen_[id];
      case ComponentType::ENVIRONMENTAL:
        return ptr_env_[id];
      case ComponentType::VERTICAL:
        return ptr_vert_[id];
      case ComponentType::TOTAL:
        return ptr_tot_[id];
    }
  }

  inline double* operator()(std::size_t id, ComponentType type) noexcept {
    switch (type) {
      case ComponentType::GENETIC:
        return ptr_gen_[id];
      case ComponentType::ENVIRONMENTAL:
        return ptr_env_[id];
      case ComponentType::VERTICAL:
        return ptr_vert_[id];
      case ComponentType::TOTAL:
        return ptr_tot_[id];
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
  std::vector<double*> ptr_gen_;
  std::vector<double*> ptr_env_;
  std::vector<double*> ptr_vert_;
  std::vector<double*> ptr_tot_;
  std::vector<bool> occupied_;
};

#endif // AMSIMCPP_PHENOBUF_H