#include <amsim/phenotype.h>
#include <amsim/genome.h>

#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <optional>
#include <algorithm>
#include <stdexcept>

#if defined(__APPLE__) && defined(USE_BLAS)
  #include <Accelerate/Accelerate.h>
#elif defined(__linux__) && defined(USE_BLAS)
  #include <cblas.h>
#endif

namespace amsim {
  PhenoBuf::PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno)
    : n_ind_(n_ind),
      n_pheno_(n_pheno) {
    // allocate and initialise buffers to storage component values
    buffer_.resize(3 * n_ind_ * n_pheno_);
    ptr_gen_.resize(n_pheno_);
    ptr_env_.resize(n_pheno_);
    ptr_vert_.resize(n_pheno_);
    occupied_.resize(n_pheno_);
    
    // setup pointers for genetic components of phenotypes
    for (std::size_t id = 0; id < n_pheno; id++) {
      ptr_gen_[id] = &buffer_[3 * id * n_ind_];
      ptr_env_[id] = &buffer_[3 * id * n_ind_ + n_ind_];
      ptr_vert_[id] = &buffer_[3 * id * n_ind_ + 2 * n_ind_];
    }
  }

  std::optional<std::size_t> PhenoBuf::unoccupied() const {
    auto it = std::find(occupied_.begin(), occupied_.end(), false);
    if (it != occupied_.end()) {
      std::size_t id = std::distance(occupied_.begin(), it);
      return id;
    }
    return std::nullopt;
  }

  void PhenoBuf::occupy(std::size_t id) {
    occupied_[id] = true;
  }

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
    if (h2_gen_ + h2_env_ + h2_vert_ != 1.0)
      throw std::runtime_error("sum of phenotype component variances must equal one");

    if (!id.has_value()) id = buf.unoccupied();
    if (!id.has_value()) throw std::runtime_error("all buffer slots occupied");
    buf.occupy(id.value());

    ptr_gen_ = buf(id.value(), ComponentType::GENETIC);
    ptr_env_ = buf(id.value(), ComponentType::ENVIRONMENTAL);
    ptr_vert_ = buf(id.value(), ComponentType::VERTICAL);
  }

  void Phenotype::score_bitwise(Genome& genome) {
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
          ptr_gen_[ind] += 1.0 * loc_effect;
        }

        for (std::uint64_t m = HOM; m; m &= (m - 1)) {
          unsigned t = static_cast<unsigned>(__builtin_ctzll(m));
          std::size_t ind = (word << 6) + t;
          ptr_gen_[ind] += 2.0 * loc_effect;
        }
      }
    }

    for (std::size_t ind = 0; ind < n_ind; ind++) 
      ptr_gen_[ind] -= global_centre;
  }
 
  #if defined(USE_BLAS)
  void Phenotype::score_tiled64(Genome& genome) {
    if (genome.H0().view() != HaploView::LOC_MAJOR)
      throw std::runtime_error("Phenotype::score_tiled64: require LOC_MAJOR view.");

    HaploBuf& H0 = genome.H0();
    HaploBuf& H1 = genome.H1();

    const std::size_t n_ind = H0.n_ind();
    const std::size_t n_words = H0.n_words();
    const std::size_t n_causal_loc = loci_.size();

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
        ptr_gen_[tile_start + k] = out[k];
    }
  }
  #endif

  void Phenotype::score(Genome& genome) {
    #ifdef USE_BLAS
      score_tiled64(genome);
    #else
      score_bitwise(genome);
    #endif
  }
}
