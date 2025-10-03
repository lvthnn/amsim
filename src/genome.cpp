#include <amsim/genome.hpp>
#include <amsim/phenotype.hpp>
#include <amsim/utils.hpp>
#include <amsim/rng.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>
#include <stdexcept>

namespace amsim {
  HapMat::HapMat(std::size_t n_ind, std::size_t n_loc)
    : n_ind_(n_ind),
      n_loc_(n_loc),
      n_rows_((n_ind + 63) & ~std::size_t(63)),
      n_words_((n_loc + 63) / 64),
      view_(HapMatView::IND_MAJOR) {
    data_.resize(n_rows_ * n_words_);
  }

  void HapMat::transpose() noexcept {
    std::size_t ct, cb, r;
    std::size_t n_tiles_row = n_rows_ / 64;
    std::size_t n_tiles_col = n_words_;
    std::uint64_t tile_cnt[64];

    for (ct = 0; ct < n_tiles_row; ct++) {
      for (cb = 0; cb < n_tiles_col; cb++) {
        for (r = 0; r < 64; r++)
          tile_cnt[r] = (*this)(ct * 64 + r, cb);
      utils::bitmatrix_transpose(tile_cnt);
        for (r = 0; r < 64; ++r)
          (*this)(ct * 64 + r, cb) = tile_cnt[r];
      }
    }

    if (view_ == HapMatView::IND_MAJOR) {
      n_rows_  = (n_loc_ + 63) & ~size_t(63);
      n_words_ = (n_ind_ + 63) / 64;
      view_    = HapMatView::LOC_MAJOR;
    } else {
      n_rows_  = (n_ind_ + 63) & ~size_t(63);
      n_words_ = (n_loc_ + 63) / 64;
      view_    = HapMatView::IND_MAJOR;
    }
  }

	Genome::Genome(std::size_t n_ind, std::size_t n_loc,
								 std::vector<double> v_mut, std::vector<double> v_rec,
								 std::vector<double> v_maf, std::uint64_t rng_seed)
		: v_mut_(std::move(v_mut)),
			v_rec_(std::move(v_rec)),
			v_maf_(std::move(v_maf)),
      bw_(rng_seed),
      H0_(n_ind, n_loc),
      H1_(n_ind, n_loc) {};

  uint64_t Genome::gam_word_(std::size_t ind, std::size_t word) noexcept {
    std::uint64_t par = bw_.sample();
    bool par0 = bw_.coinflip();
    par ^= par << 1;
    par ^= par << 2;
    par ^= par << 4;
    par ^= par << 8;
    par ^= par << 16;
    par ^= par << 32;
    if (par0) par = ~par;
    return (par & H1_(ind, word)) | (~par & H0_(ind, word));
  }

	void Genome::generate_haplotypes() noexcept {
    bw_.set_prob(0.5);

    std::size_t n_ind = H0_.n_ind();
    std::size_t n_words = H0_.n_words();

    for (size_t ind = 0; ind < n_ind; ind++) {
      std::uint64_t* word0 = H0_.rowptr(ind);
      std::uint64_t* word1 = H1_.rowptr(ind);

      for (std::size_t word = 0; word < n_words; ++word) {
        word0[word] = bw_.sample();
        word1[word] = bw_.sample();
      }
    }
	}

  void Genome::transpose() {
    H0_.transpose();
    H1_.transpose();
  }

  void Genome::compute_mafs() {
    if (H0_.view() == HapMatView::IND_MAJOR)
      throw std::runtime_error("Compute MAFs in locus-major view.");

    std::size_t n_ind = H0_.n_ind();
    std::size_t n_loc = H0_.n_loc();

		const std::size_t n_bloc_ind = (n_ind + 63) / 64;

		for (std::size_t loc = 0; loc < n_loc; ++loc) {
      std::size_t ct_loc = 0;
			for (std::size_t bloc = 0; bloc < n_bloc_ind; ++bloc) {
				if (bloc == n_bloc_ind - 1 && (n_ind % 64)) {
          std::uint64_t mask = (1ULL << (n_ind % 64)) - 1ULL;
					ct_loc += __builtin_popcountll(H0_(loc, bloc) & mask);
					ct_loc += __builtin_popcountll(H1_(loc, bloc) & mask);
				} else {
					ct_loc += __builtin_popcountll(H0_(loc, bloc));
					ct_loc += __builtin_popcountll(H1_(loc, bloc));
				}
			}
			v_lmaf_[loc] = static_cast<double>(ct_loc) / (2.0 * n_ind);
		}
  }

	void Genome::compute_stats() {
    if (H0_.view() == HapMatView::IND_MAJOR)
      throw std::runtime_error("Compute stats in loc-major view.");

    std::size_t n_ind = H0_.n_ind();
    std::size_t n_loc = H0_.n_loc();

    const std::size_t n_bloc_ind = (n_ind + 63) / 64;

    for (std::size_t loc = 0; loc < n_loc; ++loc) {
      std::size_t hom = 0;
      std::size_t het = 0;
      for (std::size_t bloc = 0; bloc < n_bloc_ind; ++bloc) {
        hom += __builtin_popcountll(H0_(loc, bloc) & H1_(loc, bloc));
        het += __builtin_popcountll(H0_(loc, bloc) ^ H1_(loc, bloc));
      }
      double mean_loc   = static_cast<double>(2.0 * hom + het) / n_ind;
      double mean_sqloc = static_cast<double>(4.0 * hom + het) / n_ind;
      double var_loc    = mean_sqloc - mean_loc * mean_loc;

      v_lmean_[loc] = mean_loc;
      v_lvar_[loc]  = var_loc;
    }
	}

  void Genome::score(Phenotype& phenotype) {
    // assume we already have processed stats
    if (H0_.view() == HapMatView::LOC_MAJOR)
      throw std::runtime_error("Score phenotypes in ind-major view.");
    
    std::size_t n_ind = H0_.n_ind();
    std::vector<double> values(n_ind);

    for (std::size_t ind = 0; ind < n_ind; ind++) {
      for (std::size_t loc = 0; loc < phenotype.loci.size(); loc++) {
        std::size_t word = n_ind / 64;
        std::size_t offset = n_ind % 64;

        double eff = phenotype.effects[loc];
        double eff_std[3];

        for (std::size_t i = 0; i < 3; i++)
          eff_std[i] = eff * (static_cast<double>(i) - v_lmean_[loc]) / v_lvar_[loc];
      }
    }

  }

	void Genome::update(std::vector<std::size_t> matching) {
		if (H0_.view() == HapMatView::LOC_MAJOR)
      throw std::runtime_error("Update in ind-major view.");

    for (std::size_t el = 0; el < matching.size(); el++)
      matching[el] -= 1;

    const std::size_t n_ind = H0_.n_ind();
    const std::size_t n_words = H0_.n_words();
    const std::size_t n_pairs = n_ind / 2;
    bw_.set_prob(0.5);

    for (std::size_t pair = 0; pair < n_pairs; pair++) {
      for (std::size_t word = 0; word < n_words; word++) {
        std::size_t fpair = pair + (n_ind / 2);
        H0_(pair, word) = gam_word_(pair, word);
        H1_(pair, word) = gam_word_(pair, word);
        H0_(fpair, word) = gam_word_(fpair, word);
        H1_(fpair, word) = gam_word_(fpair, word);
      }
    }
	}
}
