#ifndef AMSIMCPP_PHENOTYPE_H
#define AMSIMCPP_PHENOTYPE_H

#pragma once

#include <amsim/component_type.h>
#include <amsim/genome.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace amsim {

class Phenotype {
 public:
  Phenotype(
      PhenoBuf& buf,
      PhenoArch& arch,
      std::string name,
      const double h2_gen,
      const double h2_env,
      const double h2_vert,
      const std::optional<std::size_t> id = std::nullopt);

  inline const std::string name() const noexcept { return name_; }
  inline const std::vector<std::size_t> loci() const& noexcept { return loci_; }
  inline std::size_t n_ind() const noexcept { return n_ind_; }
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

  inline const double* operator()(ComponentType type) const {
    switch (type) {
      case ComponentType::GENETIC:
        return ptr_gen_;
      case ComponentType::ENVIRONMENTAL:
        return ptr_env_;
      case ComponentType::VERTICAL:
        return ptr_vert_;
      case ComponentType::TOTAL:
        return ptr_tot_;
    }
    __builtin_unreachable();
  }

  inline double comp_mean(ComponentType type) const {
    return comp_means_[type];
  }

  inline double comp_var(ComponentType type) const { return comp_vars_[type]; }

  inline void transmit_vert(std::vector<std::size_t> matching) {
    if (h2_vert_ == 0.0) return;
    const double scale = std::sqrt(h2_vert_);
    const std::size_t n_sex = n_ind_ / 2;

    for (std::size_t ind = 0; ind < n_sex; ++ind) {
      ptr_vert_[ind] = scale * (ptr_tot_[ind] + ptr_tot_[matching[ind]]);
      ptr_vert_[matching[ind]] = ptr_vert_[ind];
    }
  }

  inline void score_tot() {
    for (std::size_t ind = 0; ind < n_ind_; ++ind)
      ptr_tot_[ind] = ptr_gen_[ind] + ptr_env_[ind] + ptr_vert_[ind];
  }

  void score_bitwise(Genome& genome);
  void score_tiled(Genome& genome);
  void score(Genome& genome);
  void compute_stats();

 private:
  const std::string name_;
  const std::size_t id_;
  const std::size_t n_ind_;
  const std::vector<std::size_t> loci_;
  const std::vector<double> loc_effects_;
  const double h2_gen_;
  const double h2_env_;
  const double h2_vert_;

  double* ptr_gen_;
  double* ptr_env_;
  double* ptr_vert_;
  double* ptr_tot_;

  std::array<double, 4> comp_means_;
  std::array<double, 4> comp_vars_;
};

using PhenotypeList = std::vector<Phenotype>;

}  // namespace amsim

#endif  // AMSIMCPP_PHENOTYPE_H
