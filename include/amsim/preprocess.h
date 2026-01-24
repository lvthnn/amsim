#pragma once 

#include <amsim/params.h>

namespace amsim::preprocess {

// Preprocesses a user-supplied LD correlation matrix in order to allow
// LD-aware founder genotype generation
void calibrate_ld_matrix(Params& params);

// Optimises the causal locus assignment of phenotypes to induce a particular
// correlation structure
void optimise_phenotype_arch(Params& params);

} // namespace amsim
