#ifndef AMSIMCPP_HAPLOBUF_H
#define AMSIMCPP_HAPLOBUF_H

#include <amsim/haploview.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace amsim {

class HaploBuf {
 public:
  HaploBuf(std::size_t n_ind, std::size_t n_loc);

  std::size_t n_ind() const noexcept { return n_ind_; }
  std::size_t n_loc() const noexcept { return n_loc_; }
  std::size_t n_rows() const noexcept { return n_rows_; }
  std::size_t n_words() const noexcept { return n_words_; }
  HaploView view() const noexcept { return view_; }

  std::uint64_t& operator()(std::size_t i, std::size_t j) noexcept {
    return buf_[(i * n_words_) + j];
  }

  const std::uint64_t& operator()(std::size_t i, std::size_t j) const noexcept {
    return buf_[(i * n_words_) + j];
  }

  const std::uint64_t* rowptr(std::size_t i) const noexcept {
    return &buf_[i * n_words_];
  }

  std::uint64_t* rowptr(std::size_t i) noexcept { return &buf_[i * n_words_]; }

  void transpose() noexcept;

 private:
  const std::size_t n_ind_;
  const std::size_t n_loc_;
  std::size_t n_rows_;
  std::size_t n_words_;
  std::vector<uint64_t> buf_;
  HaploView view_;
};

}  // namespace amsim

#endif  // AMSIMCPP_HAPLOBUF_H
