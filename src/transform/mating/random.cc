#include <amsim/core/rng.h>
#include <amsim/core/utils.h>
#include <amsim/transform/mating/random.h>

#include <algorithm>
#include <numeric>

namespace amsim {

void RandomMating::randomiseMatching() {
  std::iota(match_cur_.begin(), match_cur_.end(), 0);
  for (std::size_t el = 0; el < n_sex_; ++el) {
    std::size_t le = rng::UniformIntRange::sample(el, n_sex_);
    std::swap(match_cur_[el], match_cur_[le]);
  }
}

void RandomMating::operator()(State& state) {
  randomiseMatching();
  state.matching() = match_cur_;
  for (std::size_t ind = 0; ind < n_sex_; ++ind)
    state.inv_matching()[ind] = state.matching()[ind];
}

} // namespace amsim
