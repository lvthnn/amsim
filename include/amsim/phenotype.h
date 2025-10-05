#pragma once
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include <amsim/genome.h>

namespace amsim {
  class PhenoBuffer {
  public:
    PhenoBuffer(const std::size_t n_ind, const std::size_t n_pheno,
                const std::size_t n_comps);
    inline double* ptr_gen(std::size_t id) const noexcept { return ptr_gen_[id]; }
    inline double* ptr_env(std::size_t id) const noexcept { return ptr_env_[id]; }
    inline double* ptr_vert(std::size_t id) const noexcept { return ptr_vert_[id]; }
    
  private:
    std::vector<double> buffer_;
    std::vector<double*> ptr_gen_;
    std::vector<double*> ptr_env_;
    std::vector<double*> ptr_vert_;
  };

  class Phenotype {
  public:
		Phenotype(const std::string name, const std::vector<size_t>& loci,
              const double h2);

    inline std::string name() const noexcept { return name_; }
    inline std::vector<std::size_t> loci() const& noexcept { return loci_; }
    inline double h2() const noexcept { return h2_; }

    // refactor this into a shared buffer which exposes pointers to phenotypes
    // which are stored as fields
    inline std::vector<double> values_gen() const& noexcept { return values_gen_; }
    inline std::vector<double> values_env() const& noexcept { return values_env_; }
    inline std::vector<double> values() const& noexcept { return values_; }

    void score_bitwise(Genome& genome);
    void score_tiled64(Genome& genome);
    void score(Genome& genome);
  private:
    const std::string name_;
    const std::vector<std::size_t> loci_;
    const std::vector<std::uint64_t> word_masks_;
    const double h2_;

    std::vector<double> loc_effects_;
    std::unordered_map<std::size_t, std::uint64_t> loc_mask_;

    std::vector<double> values_gen_;
    std::vector<double> values_env_;
    std::vector<double> values_;
  };
}
