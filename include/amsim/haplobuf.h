#ifndef AMSIMCPP_HAPLOBUF_H
#define AMSIMCPP_HAPLOBUF_H

#include <amsim/haploview.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace amsim {

/// @brief Class to store and manage phased haplotype data
///
/// HaploBuf implements an optimised bit-packed data structure for phased
/// haplotype data for the representation of genetic data in simulations. Two
/// views are provided to said data, in the form of an individual-major and a
/// locus-major view, respectively. The former presents the underlying data with
/// rows representing individuals and column values representing 64 bit-packed
/// locus dosages whereas the latter displays the opposite, that is, with each
/// column storing 64 individual genotypes for loci which are indicated in the
/// rows. This allows for highly optimised operations, such as the computation
/// of minor allele frequencies (MAF) and locus summary statistics.
class HaploBuf {
 public:
  /// @brief Create new HaploBuf instance
  ///
  /// Allocates a std::uint64_t buffer of size n_rows * n_words, where `n_rows`
  /// is the number of loci padded to the nearest multiple of 64, and
  /// `n_words` is the number of 64-bit blocks needed to represent all of the
  /// individual genotype dosages.
  ///
  /// @param n_ind Number of individuals to simulate in the population
  /// @param n_loc Number of genetic loci to model in the population
  HaploBuf(std::size_t n_ind, std::size_t n_loc);

  /// @brief Return the number of individuals in the population
  /// @return The number of individuals in the population
  std::size_t n_ind() const noexcept { return n_ind_; }

  /// @brief Return the number of loci modelled in the population
  /// @return The number of loci modelled in the population
  std::size_t n_loc() const noexcept { return n_loc_; }

  /// @brief Return the current number of rows in the buffer
  /// @return The current number of rows in the buffer
  std::size_t n_rows() const noexcept { return n_rows_; }

  /// @brief Return the number of words (number of columns) in the buffer
  /// @return The number of words in the buffer
  std::size_t n_words() const noexcept { return n_words_; }

  /// @brief Return the current layout of the buffer
  /// @return The current value of HaploView enumeration type
  HaploView view() const noexcept { return view_; }

  /// @brief Access the buffer data at a particular site
  /// @return The bit value at the jth column of the ith row in the current view
  std::uint64_t& operator()(std::size_t i, std::size_t j) noexcept {
    return buf_[(i * n_words_) + j];
  }

  /// @brief Access the buffer data at a particular site
  /// @param i Row index for accessing data
  /// @param j Column index for accessing data
  /// @return The (constant-access) bit value at the jth column of the ith row
  ///   in the current view
  const std::uint64_t& operator()(std::size_t i, std::size_t j) const noexcept {
    return buf_[(i * n_words_) + j];
  }

  /// @brief Access the pointer to a specified row of the buffer
  /// @param i Index of the desired row
  /// @return The pointer to the desired row to access in the buffer
  std::uint64_t* rowptr(std::size_t i) noexcept { return &buf_[i * n_words_]; }

  /// @brief Access the pointer to a particular row of the buffer
  /// @param i Index of the desired row
  /// @return The (constant-access) pointer to the desired to access row in the
  ///   buffer
  const std::uint64_t* rowptr(std::size_t i) const noexcept {
    return &buf_[i * n_words_];
  }

  /// @brief Transpose the underlying bit-matrix between views
  ///
  /// Performs block-wise transposition of the haplotype buffer, i.e., from a
  /// locus-major view to an individual-major view and vice versa.
  void transpose() noexcept;

 private:
  const std::size_t n_ind_;    ///< Number of individuals
  const std::size_t n_loc_;    ///< Number of loci
  std::size_t n_rows_;         ///< Number of rows in the buffer
  std::size_t n_words_;        ///< Number of 64-bit words per row
  std::vector<uint64_t> buf_;  ///< Bit-packed buffer storing haplotype data
  HaploView view_;             ///< Current layout view of the buffer
};

}  // namespace amsim

#endif  // AMSIMCPP_HAPLOBUF_H
