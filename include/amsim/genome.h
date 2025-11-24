#ifndef AMSIMCPP_GENOME_H
#define AMSIMCPP_GENOME_H

#include <amsim/haplobuf.h>
#include <amsim/haploview.h>
#include <amsim/rng.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace amsim {

/// @brief Class to manage and update population genotypes
///
/// Genome stores phased haplotype data in the form of two haplotype buffers,
/// which are paternal and maternal, respectively. It is also responsible for
/// generating the initial haplotypes used to start the simulation, computing
/// locus summary statistics, and updating the population genotypes to produce
/// the next generation.
///
/// @see HaploBuf
class Genome {
 public:
  /// @brief Create a new Genome instance
  ///
  /// @param n_ind Number of individuals to simulate in the population
  /// @param n_loc Number of genetic loci to model in the population
  /// @param v_mut Vector of mutation probabilities for genetic loci
  /// @param v_rec Vector of recombination probabilities for genetic loci
  /// @param v_maf vector of locus minor allele frequencies (MAFs)
  /// @param rng A seeded RNG instance
  ///
  /// @details Allocates two haplotype buffer (`HaploBuf`) instances where
  ///   haplotypes are stored in bit-packed layout initially in locus-major
  ///   layout, i.e., where each row designates a genetic locus and a column
  ///   value comprises 64 individual genotypes. Rows are padded to a multiple
  ///   of 64 to allow for efficient block-transposition of the bit matrices
  ///   used to represent haplotypes.
  ///
  /// @see HaploBuf
  Genome(
      std::size_t n_ind,
      std::size_t n_loc,
      std::vector<double> v_mut,
      std::vector<double> v_rec,
      std::vector<double> v_maf,
      const rng::Xoshiro256ss& rng);

  /// @brief Return locus means within the current generation
  /// @return A vector of locus means within the current generation
  std::vector<double>& v_lmean() noexcept { return v_lmean_; }

  /// @brief Return locus variances within the current generation
  /// @return A vector of locus variances within the current generation
  std::vector<double>& v_lvar() noexcept { return v_lvar_; }

  /// @brief Return locus MAFs within the current generation
  /// @return A vector locus MAFs within the current generation
  std::vector<double>& v_lmaf() noexcept { return v_lmaf_; }

  /// @brief Return the mean of a specified locus within the current generation
  /// @param loc The integer index of the targeted locus
  /// @return A double value representing the locus mean
  double v_lmean(std::size_t loc) const noexcept { return v_lmean_[loc]; }

  /// @brief Return the variance of a specified locus within the current
  ///   generation
  /// @param loc The integer index of the targeted locus
  /// @return A double value representing the locus variance
  double v_lvar(std::size_t loc) const noexcept { return v_lvar_[loc]; }

  /// @brief Return the variance of a specified locus within the current
  ///   generation
  /// @param loc The integer index of the targeted locus
  /// @return A double value representing the locus variance
  double v_lmaf(std::size_t loc) const noexcept { return v_lmaf_[loc]; }

  /// @brief Return the number of individuals in the population
  /// @return The number of individuals in the population
  std::size_t n_ind() const noexcept { return h0_.n_ind(); }

  /// @brief Return the number of loci modelled in the population
  /// @return The number of loci modelled in the population
  std::size_t n_loc() const noexcept { return h0_.n_loc(); }

  /// @brief Return the current layout of the haplotype buffers
  /// @return Value of the HaploView members of the haplotype buffers
  HaploView view() const noexcept { return h0_.view(); }

  /// @brief Return the paternal haplotype buffer
  /// @return HaploBuf reference pointing to paternal haplotype buffer
  HaploBuf& h0() noexcept { return h0_; }

  /// @brief Return the maternal haplotype buffer
  /// @return HaploBuf reference pointing to maternal haplotype buffer
  HaploBuf& h1() noexcept { return h1_; }

  /// @brief Generates founder haplotypes
  void generate_haplotypes() noexcept;

  /// @brief Transpose the haplotype buffers
  /// @see HaploBuf::transpose()
  void transpose() noexcept;

  /// @brief Compute locus minor allele frequncies
  void compute_mafs();

  /// @brief Compute locus summary statistics
  void compute_stats();

  /// @brief Advance the population using a mate matching
  void update(std::vector<std::size_t> matching);

 private:
  const std::vector<double> v_mut_;  ///< Mutation probabilities
  const std::vector<double> v_rec_;  ///< Recombination probabilities
  const std::vector<double> v_maf_;  ///< Minor allele frequencies
  std::vector<double> v_lmean_;      ///< Locus means
  std::vector<double> v_lvar_;       ///< Locus variances
  std::vector<double> v_lmaf_;       ///< Locus MAFs
  rng::BW16 bw_;  ///< Bernoulli word generator for sampling alleles
  HaploBuf h0_;   ///< Paternal haplotype buffer
  HaploBuf h1_;   ///< Maternal haplotype buffer

  /// @brief Generate a gamete word from two parental haplotype words
  ///
  /// @param ind_h0 Paternal haplotype word index
  /// @param ind_h1 Maternal haplotype word index
  /// @return A 64-bit word representing the gamete genotype
  uint64_t gamWord(std::uint64_t ind_h0, std::uint64_t ind_h1) noexcept;
};

}  // namespace amsim

#endif  // AMSIMCPP_GENOME_H
