#include <amsim/params.h>
#include <amsim/state.h>

namespace amsim {

State build_state(const Params& params) {
  return State{
    .genos = {genome::GenoBuf(params), genome::GenoBuf(params)},
    .phenos = {phenome::PhenoBuf(params), phenome::PhenoBuf(params)},
    .matchings = {
        std::vector<std::size_t>(params.geno.n_ind / 2),
        std::vector<std::size_t>(params.geno.n_ind / 2)}
  };
}

}
