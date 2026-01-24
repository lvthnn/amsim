#include <amsim/state.h>
#include <amsim/params.h>

#include  <amsim/initialise.h>

namespace amsim::genome {

void HaplotypeGeneratorIID::generate_haplotypes(GenoBuf& buf) {
  if (buf.view() != HaploView::LOC_MAJOR)
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

} // namespace amsim::genome
