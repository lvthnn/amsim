//------------------------------------------------------------------------------
// amsimcpp : phenotype.h
//------------------------------------------------------------------------------

#ifndef AMSIMCPP_PHENOTYPE_H
#define AMSIMCPP_PHENOTYPE_H

#pragma once
#include <cstddef>
#include <string>
#include <vector>
#include <optional>

#include <amsim/genome.h>
#include <amsim/phenobuf.h>
#include <amsim/componenttype.h>

namespace amsim {

  class Phenotype {
  public:
		Phenotype(PhenoBuf& buf, std::string name,
              std::vector<size_t>& loci, double h2_gen, double h2_vert,
              std::optional<std::size_t> id = std::nullopt);

    inline const std::string name() const noexcept { return name_; }
    inline const std::vector<std::size_t> loci() const& noexcept { return loci_; }
    inline double h2_gen() const noexcept { return h2_gen_; }
    inline double h2_vert() const noexcept { return h2_vert_; }
    inline double h2_env() const noexcept { return h2_env_; }

    inline const double& operator()(std::size_t id, ComponentType type) const {
      if (id >= n_ind_)
        throw std::runtime_error("attempting out-of-bounds access of phenotype");
      switch (type) {
        case ComponentType::GENETIC:
          return ptr_gen_[id];
        case ComponentType::ENVIRONMENTAL:
          return ptr_env_[id];
        case ComponentType::VERTICAL:
          return ptr_vert_[id];
        case ComponentType::TOTAL:
          return ptr_tot_[id];
      }
    }

    inline const double& operator()(std::size_t id) const {
      if (id >= n_ind_)
        throw std::runtime_error("attempting out-of-bounds access of phenotype");
      return ptr_tot_[id];
    }

    // @TODO: Rework this function. We need to access the mate matching of the
    //        generation and compute the sum of the two parental phenotypes
    //        weighted by some additional slop.
    inline void transmit_vert() {
      const double scale = std::sqrt(h2_vert_);
      std::copy(ptr_tot_, ptr_tot_ + n_ind_, ptr_vert_);
      for (std::size_t ind = 0; ind < n_ind_; ind++)
        ptr_vert_[ind] *= scale;
    }

    inline void score_tot() {
      for (std::size_t ind = 0; ind < n_ind_; ind++)
        ptr_tot_[ind] = ptr_gen_[ind] + ptr_env_[ind] + ptr_vert_[ind];
    }

    void score_bitwise(Genome& genome) const;
    void score_tiled64(Genome& genome) const;
    void score(Genome& genome) const;

  private:
    const std::string name_;
    const std::size_t n_ind_;
    const std::vector<std::size_t> loci_;
    const std::vector<double> loc_effects_;
    const double h2_gen_;
    const double h2_vert_;
    const double h2_env_;

    double* ptr_gen_;
    double* ptr_env_;
    double* ptr_vert_;
    double* ptr_tot_;
  };
}

#endif // AMSIMCPP_PHENOTYPE_H