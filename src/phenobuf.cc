#include <cstddef>
#include <vector>
#include <algorithm>
#include <optional>

#include <amsim/phenobuf.h>

namespace amsim {
  PhenoBuf::PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno)
    : n_ind_(n_ind),
      n_pheno_(n_pheno) {
    buf_.resize(4 * n_ind_ * n_pheno_);
    occupied_.resize(n_pheno_, false);
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
