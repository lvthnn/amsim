//
// Created by Kári Hlynsson on 7.10.2025.
//

#include <cstddef>
#include <cstdint>

#include <amsim/haplobuf.h>
#include <amsim/utils.h>

namespace amsim {

  HaploBuf::HaploBuf(std::size_t n_ind, std::size_t n_loc)
    : n_ind_(n_ind),
      n_loc_(n_loc),
      n_rows_((n_loc + 63) & ~std::size_t(63)),
      n_words_((n_ind + 63) / 64),
      view_(HaploView::LOC_MAJOR) {
    data_.resize(n_rows_ * n_words_);
  }

  void HaploBuf::transpose() noexcept {
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

    if (view_ == HaploView::IND_MAJOR) {
      n_rows_  = (n_loc_ + 63) & ~std::size_t(63);
      n_words_ = (n_ind_ + 63) / 64;
      view_    = HaploView::LOC_MAJOR;
    } else {
      n_rows_  = (n_ind_ + 63) & ~std::size_t(63);
      n_words_ = (n_loc_ + 63) / 64;
      view_    = HaploView::IND_MAJOR;
    }
  }

}