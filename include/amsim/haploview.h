#ifndef AMSIMCPP_HAPLOVIEW_H
#define AMSIMCPP_HAPLOVIEW_H

namespace amsim {

/// @brief Enumeration type representing haplotype buffer data layout
///
/// The HaploView enumeration type represents the two different layouts employed
/// by the haplotype buffer to store phased haplotype data. The former is called
/// the locus-major type, where the row indices designate loci and the values in
/// the columns comprise the genotype dosages of 64 individuals at once at the
/// corresponding locus. The latter is called the individual major view where
/// rows index individuals and the columns contain 64 genotype dosages of 64
/// loci at once.
enum class HaploView : bool {
  LOC_MAJOR,  ///< Locus-major view
  IND_MAJOR   ///< Individual-major view
};

}  // namespace amsim

#endif
