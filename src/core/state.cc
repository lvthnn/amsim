#include <amsim/core/params.h>
#include <amsim/core/state.h>

namespace amsim {

State build_state(const Params& params) {
  return State{
    .genos = {GenoBuf(params), GenoBuf(params)},
    .phenos = {PhenoBuf(params), PhenoBuf(params)},
    .matchings = {
        std::vector<std::size_t>(params.geno.n_ind / 2),
        std::vector<std::size_t>(params.geno.n_ind / 2)},
    .inv_matchings = {
        std::vector<std::size_t>(params.geno.n_ind / 2),
        std::vector<std::size_t>(params.geno.n_ind / 2)}
  };
}

}
