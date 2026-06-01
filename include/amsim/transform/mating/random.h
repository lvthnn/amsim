#pragma once

#include <amsim/core.h>
#include <amsim/data.h>
#include <amsim/transform.h>

#include <algorithm>
#include <numeric>

namespace amsim {

class RandomMating {
 public:
  explicit RandomMating(const Params& params)
    : n_sex_(params.geno.n_ind / 2),
      match_cur_(n_sex_) {}

  void operator()(State& state);

 private:
  std::size_t n_sex_;
  Matching match_cur_;

  void randomiseMatching();
};

inline void RandomMating::randomiseMatching() {
  std::iota(match_cur_.begin(), match_cur_.end(), 0);
  for (std::size_t el = 0; el < n_sex_; ++el) {
    std::size_t le = rng::UniformIntRange::sample(el, n_sex_);
    std::swap(match_cur_[el], match_cur_[le]);
  }
}

inline void RandomMating::operator()(State& state) {
  randomiseMatching();
  state.matching() = match_cur_;
  for (std::size_t ind = 0; ind < n_sex_; ++ind)
    state.inv_matching()[state.matching()[ind]] = ind;
}

} // namespace amsim
