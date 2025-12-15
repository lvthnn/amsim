#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#include <amsim/component_type.h>
#include <amsim/genome.h>
#include <amsim/haplobuf.h>
#include <amsim/logger.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/phenotype.h>
#include <amsim/stats.h>
#include <amsim/utils.h>

namespace amsim {

namespace {
std::size_t attach(PhenoBuf& buf, std::optional<std::size_t> id) {
  if (!id) id = buf.unoccupied();
  if (!id) throw std::runtime_error("all phenotype buffer slots occupied");
  return *id;
}
}  // namespace

Phenotype::Phenotype(
    PhenoBuf& buf,
    PhenoArch& arch,
    std::string name,
    const std::optional<std::size_t> id)
    : name_(std::move(name)),
      id_(attach(buf, id)),
      n_ind_(buf.n_ind()),
      loci_(arch.pheno_loc(id_)),
      loc_effects_(arch.pheno_effects(id_)),
      h2_gen_(arch.h2_gen(id_)),
      h2_env_(arch.h2_env(id_)),
      h2_vert_(arch.h2_vert(id_)),
      rvert_pat_(arch.rvert_pat(id_)),
      rvert_env_(arch.rvert_env(id_)),
      ptr_gen_(buf(id_, ComponentType::GENETIC)),
      ptr_env_(buf(id_, ComponentType::ENVIRONMENTAL)),
      ptr_vert_(buf(id_, ComponentType::VERTICAL)),
      ptr_tot_(buf(id_, ComponentType::TOTAL)),
      comp_means_(),
      comp_vars_() {
  if (h2_gen_ < 0 || h2_env_ < 0 || h2_vert_ < 0)
    throw std::runtime_error("phenotype component variances must be positive");
  if (h2_gen_ + h2_env_ + h2_vert_ != 1.0)
    throw std::runtime_error(
        "sum of phenotype component variances must equal one");

  buf.occupy(id_);
}

void Phenotype::transmit_vert(std::vector<std::size_t> matching) {
  if (h2_vert_ == 0.0) return;
  const std::size_t n_sex = n_ind_ / 2;

  for (std::size_t ind = 0; ind < n_sex; ++ind) {
    ptr_vert_[ind] = rvert_pat_ * (ptr_gen_[ind] + ptr_env_[ind]);
    ptr_vert_[ind] += (1 - rvert_pat_) * (ptr_gen_[n_sex + matching[ind]] +
                                          ptr_env_[n_sex + matching[ind]]);
    ptr_vert_[n_sex + matching[ind]] = ptr_vert_[ind];
  }

  double var_vert = stats::var(n_ind_, ptr_vert_, 1);
  LOG_INFO("var_vert: " + std::to_string(var_vert));

  if (h2_vert_ < var_vert)
    throw std::runtime_error("infeasible vertical variance noise variance");

  double var_noise = h2_vert_ - var_vert;
  LOG_INFO("var_noise: " + std::to_string(var_noise));

  if (!vert_lock_) {
    vert_var_ = std::sqrt(var_noise);
    vert_lock_ = true;
  }

  for (std::size_t ind = 0; ind < n_sex; ++ind) {
    std::array<double, 2> env_noise = rng::NormalPolar::two();
    ptr_vert_[ind] += vert_var_ * env_noise[0];
    ptr_vert_[n_sex + matching[ind]] += vert_var_ * env_noise[1];
  }
}

void Phenotype::score_bitwise(Genome& genome) {
  if (genome.view() != HaploView::LOC_MAJOR)
    throw std::runtime_error("Phenotype::score: requires loc-major view.");

  HaploBuf& h0 = genome.h0();
  HaploBuf& h1 = genome.h1();
  std::size_t n_words = h0.n_words();
  std::size_t n_ind = h0.n_ind();

  double global_centre = 0.0;

  for (std::size_t el = 0; el < loci_.size(); ++el) {
    const std::size_t loc = loci_[el];
    const double loc_sd = std::sqrt(genome.v_lvar(loc));
    const double loc_effect = loc_effects_[el] / loc_sd;
    const double loc_centre = loc_effects_[el] * genome.v_lmean(loc) / loc_sd;

    global_centre += loc_centre;

    // if the locus is monomorphic, skip it
    if (genome.v_lvar(loc) == 0) continue;

    for (std::size_t word = 0; word < n_words; ++word) {
      std::uint64_t het = h0(loc, word) ^ h1(loc, word);
      std::uint64_t hom = h0(loc, word) & h1(loc, word);
      std::size_t offset;
      std::size_t ind;

      for (std::uint64_t mask = het; mask; mask &= (mask - 1)) {
        offset = static_cast<std::size_t>(__builtin_ctzll(mask));
        ind = word * 64 + offset;
        ptr_gen_[ind] += 1.0 * loc_effect;
      }

      for (std::uint64_t mask = hom; mask; mask &= (mask - 1)) {
        offset = static_cast<std::size_t>(__builtin_ctzll(mask));
        ind = word * 64 + offset;
        ptr_gen_[ind] += 2.0 * loc_effect;
      }
    }
  }

  for (std::size_t ind = 0; ind < n_ind; ++ind) ptr_gen_[ind] -= global_centre;
}

void Phenotype::score_tiled(Genome& genome) {
  if (genome.h0().view() != HaploView::LOC_MAJOR)
    throw std::runtime_error("Phenotype::score: require LOC_MAJOR view.");

  HaploBuf& h0 = genome.h0();
  HaploBuf& h1 = genome.h1();

  const std::size_t n_ind = h0.n_ind();
  const std::size_t n_words = h0.n_words();
  const std::size_t n_causal_loc = loci_.size();

  std::vector<double> effects(n_causal_loc);
  std::vector<double> centres(n_causal_loc);
  std::vector<double> gd(64 * n_causal_loc);
  std::vector<double> out(64, 0.0);

  for (std::size_t el = 0; el < n_causal_loc; ++el) {
    std::size_t loc = loci_[el];
    double var = genome.v_lvar(loc);
    if (var == 0.0) {
      effects[el] = 0.0;
      centres[el] = 0.0;
      continue;
    }
    double sd = std::sqrt(var);
    effects[el] = loc_effects_[el] / sd;
    centres[el] = -loc_effects_[el] * genome.v_lmean(loc) / sd;
  }

  std::vector<double> ones(n_causal_loc, 1.0);

  for (std::size_t word = 0; word < n_words; ++word) {
    const std::size_t tile_start = word * 64;
    const std::size_t tile_size = std::min<std::size_t>(64, n_ind - tile_start);
    std::ranges::fill(gd, 0.0);
    std::ranges::fill(out, 0.0);

    for (std::size_t el = 0; el < n_causal_loc; ++el) {
      std::size_t loc = loci_[el];
      std::uint64_t hom = h0(loc, word) & h1(loc, word);
      std::uint64_t het = h0(loc, word) ^ h1(loc, word);

      for (std::size_t k = 0; k < tile_size; ++k) {
        int geno = (((hom >> k) & 1ULL) * 2) + ((het >> k) & 1ULL);
        gd[(k * n_causal_loc) + el] = effects[el] * geno + centres[el];
      }
    }

    cblas_dgemv(
        CblasRowMajor,
        CblasNoTrans,
        tile_size,
        n_causal_loc,
        1.0,
        gd.data(),
        n_causal_loc,
        ones.data(),
        1,
        0.0,
        out.data(),
        1);

    for (std::size_t k = 0; k < tile_size; ++k)
      ptr_gen_[tile_start + k] = out[k];
  }
}

void Phenotype::score(Genome& genome) { score_tiled(genome); }

void Phenotype::compute_stats() {
  // compute means and variances of all the components
  const std::vector<double> ones(n_ind_, 1.0);
  for (ComponentType comp = ComponentType::GENETIC;
       comp <= ComponentType::TOTAL;
       ++comp) {
    const double* ptr = (*this)(comp);
    comp_means_[comp] = stats::mean(n_ind_, ptr, 1);
    comp_vars_[comp] = stats::var(n_ind_, ptr, 1);

    if (comp == ComponentType::TOTAL) break;
  }
}
}  // namespace amsim
