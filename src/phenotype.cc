//------------------------------------------------------------------------------
// amsimcpp : phenotype.cc
//------------------------------------------------------------------------------

#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <optional>
#include <stdexcept>

#include <amsim/componenttype.h>
#include <amsim/haplobuf.h>
#include <amsim/genome.h>
#include <amsim/phenobuf.h>
#include <amsim/phenotype.h>

#if defined(__APPLE__) && defined(USE_BLAS)
  #include <Accelerate/Accelerate.h>
#elif defined(__linux__) && defined(USE_BLAS)
  #include <cblas.h>
#endif

namespace amsim {
  Phenotype::Phenotype(PhenoBuf& buf, std::string name,
                       std::vector<std::size_t>& loci, double h2_gen,
                       double h2_vert, std::optional<std::size_t> id)
    : name_(std::move(name)),
      n_ind_(buf.n_ind()),
      loci_(std::move(loci)),
      loc_effects_(loci_.size(), std::sqrt(h2_gen / loci_.size())),
      h2_gen_(h2_gen),
      h2_vert_(h2_vert),
      h2_env_(1 - h2_gen - h2_vert) {
    if (h2_gen_ < 0 || h2_env_ < 0 || h2_vert_ < 0)
      throw std::runtime_error("phenotype component variances must be positive");
    if (h2_gen_ + h2_env_ + h2_vert_ != 1.0)
      throw std::runtime_error("sum of phenotype component variances must equal one");

    if (!id.has_value()) id = buf.unoccupied();
    if (!id.has_value()) throw std::runtime_error("all buffer slots occupied");
    const std::size_t id_val = id.value();

    ptr_gen_  = buf(id_val, ComponentType::GENETIC);
    ptr_env_  = buf(id_val, ComponentType::ENVIRONMENTAL);
    ptr_vert_ = buf(id_val, ComponentType::VERTICAL);
    ptr_tot_  = buf(id_val, ComponentType::TOTAL);
    buf.occupy(id_val);
  }

  void Phenotype::score_bitwise(Genome& genome) const {
    if (genome.view() != HaploView::LOC_MAJOR)
      throw std::runtime_error("Phenotype::score: requires loc-major view.");

    HaploBuf& H0 = genome.H0();
    HaploBuf& H1 = genome.H1();
    std::size_t n_words = H0.n_words();
    std::size_t n_ind = H0.n_ind();

    double global_centre = 0.0;
  
    for (std::size_t el = 0; el < loci_.size(); el++) {
      const std::size_t loc = loci_[el];
      const double loc_sd = std::sqrt(genome.v_lvar(loc));
      const double loc_effect = loc_effects_[el] / loc_sd;
      const double loc_centre = loc_effects_[el] * genome.v_lmean(loc) / loc_sd;
      
      global_centre += loc_centre;

      // if the locus is monomorphic, skip it
      if (genome.v_lvar(loc) == 0) continue;
      
      for (std::size_t word = 0; word < n_words; word++) {
        std::uint64_t HET = H0(loc, word) ^ H1(loc, word);
        std::uint64_t HOM = H0(loc, word) & H1(loc, word);
        std::size_t offset, ind;

        for (std::uint64_t mask = HET; mask; mask &= (mask - 1)) {
          offset = static_cast<std::size_t>(__builtin_ctzll(mask));
          ind = word * 64 + offset;
          ptr_gen_[ind] += 1.0 * loc_effect;
        }

        for (std::uint64_t mask = HOM; mask; mask &= (mask - 1)) {
          offset = static_cast<std::size_t>(__builtin_ctzll(mask));
          ind = word * 64 + offset;
          ptr_gen_[ind] += 2.0 * loc_effect;
        }
      }
    }

    for (std::size_t ind = 0; ind < n_ind; ind++) 
      ptr_gen_[ind] -= global_centre;
  }
 
  #if defined(USE_BLAS)
  void Phenotype::score_tiled64(Genome& genome) const {
    if (genome.H0().view() != HaploView::LOC_MAJOR)
      throw std::runtime_error("Phenotype::score_tiled64: require LOC_MAJOR view.");

    HaploBuf& H0 = genome.H0();
    HaploBuf& H1 = genome.H1();

    const std::size_t n_ind = H0.n_ind();
    const std::size_t n_words = H0.n_words();
    const std::size_t n_causal_loc = loci_.size();

    std::vector<double> effects(n_causal_loc);
    std::vector<double> centres(n_causal_loc);
    std::vector<double> Gd(64 * n_causal_loc);
    std::vector<double> out(64, 0.0);

    for (std::size_t el = 0; el < n_causal_loc; ++el) {
      std::size_t loc = loci_[el];
      double var = genome.v_lvar(loc);
      if (var == 0.0) { effects[el] = 0.0; centres[el] = 0.0; continue; }
      double sd = std::sqrt(var);
      effects[el] = loc_effects_[el] / sd;
      centres[el] = -loc_effects_[el] * genome.v_lmean(loc) / sd;
    }

    std::vector<double> ones(n_causal_loc, 1.0);

    for (std::size_t word = 0; word < n_words; ++word) {
      const std::size_t tile_start = word * 64;
      const std::size_t tile_size = std::min<std::size_t>(64, n_ind - tile_start);
      std::fill(Gd.begin(), Gd.end(), 0.0);
      std::fill(out.begin(), out.end(), 0.0);

      for (std::size_t el = 0; el < n_causal_loc; ++el) {
        std::size_t loc = loci_[el];
        std::uint64_t HOM = H0(loc, word) & H1(loc, word);
        std::uint64_t HET = H0(loc, word) ^ H1(loc, word);

        for (std::size_t k = 0; k < tile_size; ++k) {
          int geno = ((HOM >> k) & 1ull) * 2 + ((HET >> k) & 1ull);
          Gd[k * n_causal_loc + el] = effects[el] * geno + centres[el];
        }
      }

      cblas_dgemv(CblasRowMajor, CblasNoTrans,
                  tile_size, n_causal_loc,
                  1.0,
                  Gd.data(), n_causal_loc,
                  ones.data(), 1,
                  0.0,
                  out.data(), 1);

      for (std::size_t k = 0; k < tile_size; ++k)
        ptr_gen_[tile_start + k] = out[k];
    }
  }
  #endif

  void Phenotype::score(Genome& genome) const {
    #ifdef USE_BLAS
      score_tiled64(genome);
    #else
      score_bitwise(genome);
    #endif
  }
}
