#include <amsim/haplobuf.h>
#include <amsim/utils.h>

#include <cstddef>
#include <cstdint>

namespace amsim {

HaploBuf::HaploBuf(std::size_t n_ind, std::size_t n_loc)
    : n_ind_(n_ind),
      n_loc_(n_loc),
      n_rows_((n_loc + 63) & ~std::size_t(63)),
      n_words_((n_ind + 63) / 64),
      view_(HaploView::LOC_MAJOR) {
  buf_.resize(n_rows_ * n_words_);
}

void HaploBuf::transpose() noexcept {
  // Source geometry (pre-transpose)
  const std::size_t src_rows = n_rows_;
  const std::size_t src_words = n_words_;

  // Destination geometry (post-transpose)
  const std::size_t dst_rows = (view_ == HaploView::LOC_MAJOR)
                                   ? ((n_ind_ + 63) & ~std::size_t(63))
                                   : ((n_loc_ + 63) & ~std::size_t(63));

  const std::size_t dst_cols_w = (view_ == HaploView::LOC_MAJOR)
                                     ? ((n_loc_ + 63) / 64)
                                     : ((n_ind_ + 63) / 64);

  // Allocate destination buffer
  std::vector<std::uint64_t> out(dst_rows * dst_cols_w, 0);

  // Tile counts: src_rows is multiple of 64 by construction
  const std::size_t src_tile_rows = src_rows / 64;
  const std::size_t src_tile_cols = src_words;

  std::uint64_t tile[64];

  for (std::size_t ct = 0; ct < src_tile_rows; ++ct) {
    for (std::size_t cb = 0; cb < src_tile_cols; ++cb) {
      for (std::size_t r = 0; r < 64; ++r) tile[r] = (*this)(ct * 64 + r, cb);

      utils::bitmatrix_transpose(tile);

      for (std::size_t r = 0; r < 64; ++r) {
        const std::size_t dst_row = cb * 64 + r;  // 64 rows per source word
        const std::size_t dst_col = ct;  // word index becomes column index
        out[dst_row * dst_cols_w + dst_col] = tile[r];
      }
    }
  }

  buf_.swap(out);
  n_rows_ = dst_rows;
  n_words_ = dst_cols_w;
  view_ = (view_ == HaploView::LOC_MAJOR) ? HaploView::IND_MAJOR
                                          : HaploView::LOC_MAJOR;
}

}  // namespace amsim
