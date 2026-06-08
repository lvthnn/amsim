// This file is part of amsim, copyright (C) 2025-2026 Kári Hlynsson.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
