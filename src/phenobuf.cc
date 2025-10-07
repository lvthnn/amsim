//------------------------------------------------------------------------------
// amsimcpp : phenobuf.cc
//------------------------------------------------------------------------------

#include <cstddef>
#include <optional>
#include <algorithm>

namespace amsim {
  PhenoBuf::PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno)
    : n_ind_(n_ind),
      n_pheno_(n_pheno) {
    buffer_.resize(4 * n_ind_ * n_pheno_);
    ptr_gen_.resize(n_pheno_);
    ptr_env_.resize(n_pheno_);
    ptr_vert_.resize(n_pheno_);
    ptr_tot_.resize(n_pheno_);
    occupied_.resize(n_pheno_);

    // setup pointers for genetic components of phenotypes
    for (std::size_t id = 0; id < n_pheno; id++) {
      ptr_gen_[id]  = &buffer_[id * n_ind_];
      ptr_env_[id]  = &buffer_[(n_pheno_ + id) * n_ind_];
      ptr_vert_[id] = &buffer_[(2 * n_pheno_ + id) * n_ind_];
      ptr_tot_[id]  = &buffer_[(3 * n_pheno_ + id) * n_ind_];
    }
  }

  std::optional<std::size_t> PhenoBuf::unoccupied() const {
    auto it = std::find(occupied_.begin(), occupied_.end(), false);
    if (it != occupied_.end()) {
      std::size_t id = std::distance(occupied_.begin(), it);
      return id;
    }
    return std::nullopt;
  }

  void PhenoBuf::occupy(std::size_t id) {
    occupied_[id] = true;
  }
}
