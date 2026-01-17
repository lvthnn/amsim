#include <amsim/component_type.h>
#include <amsim/genome.h>
#include <amsim/haploview.h>
#include <amsim/logger.h>
#include <amsim/rng.h>
#include <amsim/stats.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#if __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

namespace amsim {
Genome::Genome(
    std::size_t n_ind,
    std::size_t n_loc,
    std::vector<double> v_mut,
    std::vector<double> v_rec,
    std::vector<double> v_maf)
    : v_mut_(std::move(v_mut)),
      v_rec_(std::move(v_rec)),
      v_maf_(std::move(v_maf)),
      v_lmean_(n_loc),
      v_lvar_(n_loc),
      v_lmaf_(n_loc),
      bw_(),
      h0_(n_ind, n_loc),
      h1_(n_ind, n_loc) {};

std::array<std::uint64_t, 2> Genome::gamWord(
    std::uint64_t ind_h0,
    std::uint64_t ind_h1,
    const double* v_rec,
    const double* v_mut) noexcept {
  // Set recombination probabilities for loci in word
  bw_.set_probs(v_rec);

  // Sample a 0-1 recombination mask
  std::uint64_t par = bw_.sample();

  // Select the initial parental strand uniformly
  bool par0 = bw_.coinflip();

  // Hallis-Steele shift cumulative sum mod 2
  par ^= par << 1;
  par ^= par << 2;
  par ^= par << 4;
  par ^= par << 8;
  par ^= par << 16;
  par ^= par << 32;
  if (par0) par = ~par;

  // Set mutation probabilities for the loci
  bw_.set_probs(v_mut);
  std::uint64_t mut = bw_.sample();

  return {par, ((par & ind_h0) | (~par & ind_h1)) ^ mut};
}

void Genome::generate_haplotypes() noexcept {
  std::size_t n_loc = h0_.n_loc();
  std::size_t n_words = h0_.n_words();

  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    // Use a single probability with BW generator
    bw_.set_prob(v_maf_[loc]);

    // Row pointers for easy access
    std::uint64_t* word0 = h0_.rowptr(loc);
    std::uint64_t* word1 = h1_.rowptr(loc);

    for (std::size_t word = 0; word < n_words; ++word) {
      word0[word] = bw_.sample();
      word1[word] = bw_.sample();
    }
  }
}

void Genome::transpose() noexcept {
  h0_.transpose();
  h1_.transpose();
}

void Genome::compute_mafs() {
  if (h0_.view() == HaploView::IND_MAJOR)
    throw std::runtime_error("Compute MAFs in locus-major view.");

  std::size_t n_ind = h0_.n_ind();
  std::size_t n_loc = h0_.n_loc();

  const std::size_t n_bloc_ind = (n_ind + 63) / 64;

  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    // Total locus dosage
    std::size_t ct_loc = 0;

    for (std::size_t bloc = 0; bloc < n_bloc_ind; ++bloc) {
      if (bloc == n_bloc_ind - 1 && (n_ind % 64)) {
        // Mask to trim off padding in last block
        std::uint64_t mask = (1ULL << (n_ind % 64)) - 1ULL;
        ct_loc += __builtin_popcountll(h0_(loc, bloc) & mask);
        ct_loc += __builtin_popcountll(h1_(loc, bloc) & mask);
      } else {
        // Popcount words for fast total haplotype dosage
        ct_loc += __builtin_popcountll(h0_(loc, bloc));
        ct_loc += __builtin_popcountll(h1_(loc, bloc));
      }
    }

    // Compute the MAF
    v_lmaf_[loc] =
        static_cast<double>(ct_loc) / (2.0 * static_cast<double>(n_ind));
  }
}

void Genome::compute_stats() {
  if (h0_.view() == HaploView::IND_MAJOR)
    throw std::runtime_error("Compute stats in loc-major view.");

  std::size_t n_ind = h0_.n_ind();
  std::size_t n_loc = h0_.n_loc();

  const std::size_t n_bloc_ind = (n_ind + 63) / 64;

  for (std::size_t loc = 0; loc < n_loc; ++loc) {
    std::size_t hom = 0;
    std::size_t het = 0;
    for (std::size_t bloc = 0; bloc < n_bloc_ind; ++bloc) {
      if (bloc == n_bloc_ind - 1 && (n_ind % 64)) {
        std::uint64_t mask = (1ULL << (n_ind % 64)) - 1ULL;
        hom += __builtin_popcountll((h0_(loc, bloc) & h1_(loc, bloc)) & mask);
        het += __builtin_popcountll((h0_(loc, bloc) ^ h1_(loc, bloc)) & mask);
      } else {
        hom += __builtin_popcountll(h0_(loc, bloc) & h1_(loc, bloc));
        het += __builtin_popcountll(h0_(loc, bloc) ^ h1_(loc, bloc));
      }
    }
    double mean_loc = static_cast<double>((2.0 * hom) + het) / n_ind;
    double mean_sqloc = static_cast<double>((4.0 * hom) + het) / n_ind;
    double var_loc = mean_sqloc - (mean_loc * mean_loc);

    v_lmean_[loc] = mean_loc;
    v_lvar_[loc] = var_loc;
  }
}

void Genome::update(
    std::vector<std::size_t> matching, PhenoArch& arch, PhenoBuf& buf) {
  if (h0_.view() == HaploView::LOC_MAJOR)
    throw std::runtime_error("Update in ind-major view.");

  constexpr std::size_t INC_WORD = 64;
  const std::size_t n_ind = h0_.n_ind();
  const std::size_t n_sex = n_ind / 2;
  const std::size_t n_words = h0_.n_words();

  std::vector<std::uint64_t> nt_mm(n_words * n_sex);
  std::vector<std::uint64_t> nt_mf(n_words * n_sex);
  std::vector<std::uint64_t> nt_fm(n_words * n_sex);
  std::vector<std::uint64_t> nt_ff(n_words * n_sex);

  for (std::size_t pair = 0; pair < n_sex; ++pair) {
    std::size_t fpair = matching[pair] + n_sex;
    const double* rec_ptr = v_rec_.data();
    const double* mut_ptr = v_mut_.data();
    for (std::size_t word = 0; word < n_words; ++word) {
      std::size_t idx = (pair * n_words) + word;

      std::uint64_t male_h0 = h0_(pair, word);
      std::uint64_t male_h1 = h1_(pair, word);
      std::uint64_t female_h0 = h0_(fpair, word);
      std::uint64_t female_h1 = h1_(fpair, word);

      // male child
      auto [par_mm, trans_mm] = gamWord(male_h0, male_h1, rec_ptr, mut_ptr);
      auto [par_fm, trans_fm] = gamWord(female_h0, female_h1, rec_ptr, mut_ptr);
      h0_(pair, word) = trans_mm;  // paternal transmission to male child
      h1_(pair, word) = trans_fm;  // maternal transmission to male child

      // female child
      auto [par_mf, trans_mf] = gamWord(male_h0, male_h1, rec_ptr, mut_ptr);
      auto [par_ff, trans_ff] = gamWord(female_h0, female_h1, rec_ptr, mut_ptr);
      h0_(fpair, word) = trans_mf;  // paternal transmission to female child
      h1_(fpair, word) = trans_ff;  // maternal transmission to female child

      // non-transmitted paternal grandfather to father
      nt_mm[idx] = (~par_mm & male_h0) | (par_mm & male_h1);

      // non-transmitted maternal grandfather to mother
      nt_mf[idx] = (~par_mf & male_h0) | (par_mf & male_h1);

      // non-transmitted paternal grandmother to father
      nt_fm[idx] = (~par_fm & female_h0) | (par_fm & female_h1);

      // non-transmitted maternal grandmother to father
      nt_ff[idx] = (~par_ff & female_h0) | (par_ff & female_h1);

      rec_ptr += INC_WORD;
      mut_ptr += INC_WORD;
    }
  }

  for (std::size_t pheno = 0; pheno < arch.n_pheno(); ++pheno) {
    if (arch.h2_vert(pheno) == 0.0) continue;
    const std::vector<std::size_t>& loci = arch.pheno_loc(pheno);
    const double* effects = arch.pheno_effects(pheno);

    double* ptr_vert = buf(pheno, ComponentType::VERTICAL);
    double* ptr_env = buf(pheno, ComponentType::ENVIRONMENTAL);
    double rvert_pat = arch.rvert_pat(pheno);
    double rvert_env = arch.rvert_env(pheno);
    double rvert_mat = 1 - rvert_pat;

    std::ranges::fill_n(ptr_vert, n_ind, 0.0);

    for (std::size_t el = 0; el < loci.size(); ++el) {
      const std::size_t loc = loci[el];
      const std::size_t word = loc / 64;
      const std::size_t offset = loc % 64;

      const double loc_sd = std::sqrt(v_lvar_[loc]);
      const double loc_effect = effects[el] / loc_sd;
      const double loc_centre = loc_effect * v_lmean_[el];

      for (std::size_t ind = 0; ind < n_sex; ++ind) {
        std::size_t idx = (ind * n_words) + word;

        int nt_pat_m = ((nt_mm[idx] & (1ULL << offset)) >> offset);
        int nt_mat_m = ((nt_fm[idx] & (1ULL << offset)) >> offset);
        int nt_pat_f = ((nt_mf[idx] & (1ULL << offset)) >> offset);
        int nt_mat_f = ((nt_ff[idx] & (1ULL << offset)) >> offset);
        double snt_pat_m = (loc_effect * nt_pat_m) - loc_centre;
        double snt_mat_m = (loc_effect * nt_mat_m) - loc_centre;
        double snt_pat_f = (loc_effect * nt_pat_f) - loc_centre;
        double snt_mat_f = (loc_effect * nt_mat_f) - loc_centre;

        double nt_male = (rvert_pat * snt_pat_m) + (rvert_mat * snt_mat_m);
        double nt_female = (rvert_pat * snt_pat_f) + (rvert_mat * snt_mat_f);

        ptr_vert[ind] += (1 - rvert_env) * nt_male;
        ptr_vert[n_sex + matching[ind]] += (1 - rvert_env) * nt_female;
      }
    }

    for (std::size_t ind = 0; ind < n_sex; ++ind) {
      double par_env =
          rvert_env * ((rvert_pat * ptr_env[ind]) +
                       ((1 - rvert_pat) * ptr_env[n_sex + matching[ind]]));
      ptr_vert[ind] += par_env;
      ptr_vert[n_sex + matching[ind]] += par_env;
    }

    if (!arch.vert_lock(pheno)) {
      double var_vert = stats::var(n_ind, ptr_vert, 1);
      arch.calibrate_scale(pheno, var_vert);
    }

    cblas_dscal(n_ind, arch.vert_scale(pheno), ptr_vert, 1);
    double var_vert = stats::var(n_ind, ptr_vert, 1);
    LOG_DEBUG("var_vert: " + std::to_string(var_vert));
  }
}
}  // namespace amsim
