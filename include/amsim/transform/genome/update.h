#pragma once

#include <amsim/core/params.h>
#include <amsim/core/state.h>

#include <vector>

namespace amsim {

class UpdateGenome {
 public:
  explicit UpdateGenome(const Params& params)
      : n_ind_(params.geno.n_ind),
        n_sex_(n_ind_ / 2),
        v_rec_(params.geno.v_rec),
        v_mut_(params.geno.v_mut),
        bw_() {};
  void operator()(State& state);

 private:
  const std::size_t n_ind_;
  const std::size_t n_sex_;
  const Eigen::VectorXd& v_rec_;
  const Eigen::VectorXd& v_mut_;
  const double* ptr_rec_;
  const double* ptr_mut_;
  std::vector<std::uint64_t> transmit_chunk_;
  rng::BernoulliWord<16> bw_;

  std::array<std::uint64_t, 2> gamWord(std::uint64_t h0, std::uint64_t h1);

  void updateGenome(State& state);
  void updateNurture(State& state);
};

} // namespace amsim
