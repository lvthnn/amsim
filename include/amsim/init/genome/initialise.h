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

#include <amsim/core/params.h>
#include <amsim/core/state.h>

namespace amsim {

class HaplotypeGenerator {
 public:
  virtual ~HaplotypeGenerator() = default;
  virtual void generateHaplotypes(GenoBuf& buf) = 0;
};

// Class to initialise founder population genotypes from unlinked loci
class HaplotypeGeneratorIID : public HaplotypeGenerator {
 public:
  explicit HaplotypeGeneratorIID(const Params& params)
      : v_maf_(params.geno.v_maf), bw_() {};

  void operator()(State& state);

 private:
  const Eigen::VectorXd& v_maf_;
  rng::BernoulliWord<16> bw_;

  void generateHaplotypes(GenoBuf& buf) override;
};

inline void HaplotypeGeneratorIID::generateHaplotypes(GenoBuf& buf) {
  if (buf.view() != HaploView::LocusMajor)
    throw std::runtime_error("generate haplotypes in loc-major view");

  HaploBuf& h0 = buf.h0();
  HaploBuf& h1 = buf.h1();
  std::size_t n_loc = buf.n_loc();
  std::size_t n_words = buf.n_words();

  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    // set probability of bernoulli generator
    bw_.set_prob(v_maf_(loc));

    // row pointers for easy access
    std::uint64_t* word0 = h0.rowptr(loc);
    std::uint64_t* word1 = h1.rowptr(loc);

    for (std::size_t word = 0; word < n_words; ++word) {
      word0[word] = bw_.sample();
      word1[word] = bw_.sample();
    }
  }
}

inline void HaplotypeGeneratorIID::operator()(State& state) {
  generateHaplotypes(state.geno());
}

}  // namespace amsim
