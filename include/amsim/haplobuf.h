//------------------------------------------------------------------------------
// amsimcpp : haplobuf.h
//------------------------------------------------------------------------------

#ifndef AMSIMCPP_HAPLOBUF_H
#define AMSIMCPP_HAPLOBUF_H

#include <cstddef>
#include <cstdint>
#include <vector>

namespace amsim {
  enum class HaploView : bool {
    IND_MAJOR,
    LOC_MAJOR
  };

  class HaploBuf {
  public:
    HaploBuf(std::size_t n_ind, std::size_t n_loc);

    inline std::size_t n_ind() const noexcept { return n_ind_; }
    inline std::size_t n_loc() const noexcept { return n_loc_; }
    inline std::size_t n_rows() const noexcept { return n_rows_; }
    inline std::size_t n_words() const noexcept { return n_words_; }
    inline HaploView view() const noexcept { return view_; }

    inline std::uint64_t& operator()(std::size_t i, std::size_t j) noexcept {
      return data_[i * n_words_ + j];
    }

    inline const std::uint64_t& operator()(std::size_t i, std::size_t j) const noexcept {
      return data_[i * n_words_ + j];
    }

    inline const std::uint64_t* rowptr(std::size_t i) const noexcept {
      return &data_[i * n_words_];
    }

    inline std::uint64_t* rowptr(std::size_t i) noexcept {
      return &data_[i * n_words_];
    }

    void transpose() noexcept;

  private:
    const std::size_t n_ind_;
    const std::size_t n_loc_;
    std::size_t n_rows_;
    std::size_t n_words_;
    std::vector<uint64_t> data_;
    HaploView view_;
  };
}

#endif //AMSIMCPP_HAPLOBUF_H