#pragma once

#include <amsim/state.h>
#include <amsim/params.h>

namespace amsim {

namespace genome {

class HaplotypeGenerator {
 public:
  virtual ~HaplotypeGenerator() = default;
  virtual void generate_haplotypes(GenoBuf& buf) = 0;
};

// Class to initialise founder population genotypes from unlinked loci
class HaplotypeGeneratorIID : public HaplotypeGenerator {
 public:
  explicit HaplotypeGeneratorIID(const Params& params)
    : v_maf_(params.geno.v_maf),
      bw_() {};
  void generate_haplotypes(GenoBuf& buf) override;

 private:
  const Eigen::VectorXd& v_maf_;
  rng::BernoulliWord<16> bw_;
};


// Class to initialise founder population genotypes with linkage structure
class HaploGeneratorLD : public HaplotypeGenerator {
 public:
  explicit HaploGeneratorLD(const Params& params);
  void generate_haplotypes(GenoBuf& buf) override;

 private:
  std::vector<double> v_maf_;
  std::vector<double> v_rec_;
  std::vector<Eigen::MatrixXd> ld_blocks_;
  std::vector<std::size_t> n_loc_block_;
  std::size_t n_ld_blocks_;
};

} // namespace genome


namespace phenome {

// assign phenotype effect vectors to produce
class PhenotypeEffects {
 public:

 private:
  Eigen::VectorXd h2_gen_;
  Eigen::MatrixXd ld_cor_;
  Eigen::MatrixXd gen_cor_;
};


} // namespace phenome


} // namespace amsim
