#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>
#include <amsim/data/mating.h>

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


} // namespace amsim
