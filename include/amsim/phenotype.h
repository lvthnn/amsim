#pragma once
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include <amsim/genome.h>

namespace amsim {
  class Phenotype {
  public:
		Phenotype(const std::string name, const std::vector<size_t>& loci,
              const double h2);

    inline std::string name() const noexcept { return name_; }
    inline std::vector<std::size_t> loci() const& noexcept { return loci_; }
    inline double h2() const noexcept { return h2_; }

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
