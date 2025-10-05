#include <amsim/phenotype.h>
#include <amsim/genome.h>

#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <unordered_map>
#include <stdexcept>

#if defined(__APPLE__) && defined(USE_BLAS)
  #include <Accelerate/Accelerate.h>
#elif defined(__linux__) && defined(USE_BLAS)
  #include <cblas.h>
#endif

namespace amsim {
  
  Phenotype::Phenotype(const std::string name,
                       const std::vector<std::size_t>& loci,
                       const double h2)
    : name_(name),
      loci_(std::move(loci)),
      h2_(h2) {
    for (std::size_t el = 0; el < loci_.size(); el++) {
      std::size_t loc = loci_[el];
      std::size_t block = loc / 64;
      std::size_t offset = loc % 64;
      loc_mask_[block] |= std::uint64_t(1) << offset;
    }

    loc_effects_.resize(loci.size(), std::sqrt(h2 / loci.size()));
  }

  void Phenotype::score_bitwise(Genome& genome) {
    if (genome.view() != HapMatView::LOC_MAJOR)
      throw std::runtime_error("Phenotype::score: requires loc-major view.");

    HapMat& H0 = genome.H0();
    HapMat& H1 = genome.H1();
    std::size_t n_words = H0.n_words();
    std::size_t n_ind = H0.n_ind();

    double global_centre = 0.0;

    values_gen_.resize(n_ind, 0.0);
  
    for (std::size_t el = 0; el < loci_.size(); el++) {
      const std::size_t loc = loci_[el];
      const double loc_sd = std::sqrt(genome.v_lvar(loc));
      const double loc_effect = loc_effects_[loc] / loc_sd;
      const double loc_centre = loc_effects_[loc] * genome.v_lmean(loc) / loc_sd;
      
      global_centre += loc_centre;

      // if the locus is monomorphic, skip it
      if (genome.v_lvar(loc) == 0) continue;
      
      for (std::size_t word = 0; word < n_words; word++) {
        std::uint64_t HET = H0(loc, word) ^ H1(loc, word);
        std::uint64_t HOM = H0(loc, word) & H1(loc, word);

        for (std::uint64_t m = HET; m; m &= (m - 1)) {
          unsigned t = static_cast<unsigned>(__builtin_ctzll(m));
          std::size_t ind = (word << 6) + t;
          values_gen_[ind] += 1.0 * loc_effect;
        }

        for (std::uint64_t m = HOM; m; m &= (m - 1)) {
          unsigned t = static_cast<unsigned>(__builtin_ctzll(m));
          std::size_t ind = (word << 6) + t;
          values_gen_[ind] += 2.0 * loc_effect;
        }
      }
    }

    for (std::size_t ind = 0; ind < n_ind; ind++) 
      values_gen_[ind] -= global_centre;
  }
  
  void Phenotype::score_tiled64(Genome& genome) {
      if (genome.H0().view() != HapMatView::LOC_MAJOR)
          throw std::runtime_error("Phenotype::score_tiled64: require LOC_MAJOR view.");

      HapMat& H0 = genome.H0();
      HapMat& H1 = genome.H1();

      const std::size_t n_ind = H0.n_ind();
      const std::size_t n_words = H0.n_words();
      const std::size_t n_causal_loc = loci_.size();

      values_gen_.assign(n_ind, 0.0);

      std::vector<double> effects(n_causal_loc);
      std::vector<double> centres(n_causal_loc);

      for (std::size_t el = 0; el < n_causal_loc; ++el) {
          std::size_t loc = loci_[el];
          double var = genome.v_lvar(loc);
          if (var == 0.0) { effects[el] = 0.0; centres[el] = 0.0; continue; }
          double sd = std::sqrt(var);
          effects[el] = loc_effects_[el] / sd;
          centres[el] = -loc_effects_[el] * genome.v_lmean(loc) / sd;
      }

      std::vector<double> ones(n_causal_loc, 1.0);

      for (std::size_t b = 0; b < n_words; ++b) {
          const std::size_t tile_start = b * 64;
          const std::size_t tile_size = std::min<std::size_t>(64, n_ind - tile_start);

          std::vector<double> Gd(tile_size * n_causal_loc);

          for (std::size_t el = 0; el < n_causal_loc; ++el) {
              std::size_t loc = loci_[el];
              std::uint64_t HOM = H0(loc, b) & H1(loc, b);
              std::uint64_t HET = H0(loc, b) ^ H1(loc, b);

              for (std::size_t k = 0; k < tile_size; ++k) {
                  int geno = ((HOM >> k) & 1ull) * 2 + ((HET >> k) & 1ull);
                  Gd[k * n_causal_loc + el] = effects[el] * geno + centres[el];
              }
          }

          std::vector<double> out(tile_size, 0.0);
          cblas_dgemv(CblasRowMajor, CblasNoTrans,
                      tile_size, n_causal_loc,
                      1.0,
                      Gd.data(), n_causal_loc,
                      ones.data(), 1,
                      0.0,
                      out.data(), 1);

          for (std::size_t k = 0; k < tile_size; ++k)
              values_gen_[tile_start + k] = out[k];
      }
  }

  void Phenotype::score(Genome& genome) {
    #ifdef USE_BLAS
      score_tiled64(genome);
    #else
      score_bitwise(genome);
    #endif
  }
}
