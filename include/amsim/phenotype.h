#pragma once
#include <cstddef>
#include <string>
#include <vector>
#include <optional>

#include <amsim/genome.h>

namespace amsim {
  enum class ComponentType {
    GENETIC,
    ENVIRONMENTAL,
    VERTICAL
  };

  class PhenoBuf {
  public:
    PhenoBuf(const std::size_t n_ind, const std::size_t n_pheno);

    inline const double* operator()(std::size_t id, ComponentType type) const noexcept {
      if (type == ComponentType::GENETIC)
        return ptr_gen_[id];
      else if (type == ComponentType::ENVIRONMENTAL)
        return ptr_env_[id];
      return ptr_vert_[id];
    }

    inline double* operator()(std::size_t id, ComponentType type) noexcept {
      if (type == ComponentType::GENETIC)
        return ptr_gen_[id];
      else if (type == ComponentType::ENVIRONMENTAL)
        return ptr_env_[id];
      return ptr_vert_[id];
    }

    inline std::size_t n_ind() const noexcept { return n_ind_; }
    inline bool occupied(std::size_t id) const noexcept { return occupied_[id]; }
  
    std::optional<std::size_t> unoccupied() const;
    void occupy(std::size_t);
    
  private:
    const std::size_t n_ind_;
    const std::size_t n_pheno_;
    std::vector<double> buffer_;
    std::vector<double*> ptr_gen_;
    std::vector<double*> ptr_env_;
    std::vector<double*> ptr_vert_;
    std::vector<bool> occupied_;
  };

  class Phenotype {
  public:
		Phenotype(PhenoBuf& buf, std::string name,
              std::vector<size_t>& loci, double h2_gen, double h2_vert,
              std::optional<std::size_t> id);

    inline const std::string name() const noexcept { return name_; }
    inline const std::vector<std::size_t> loci() const& noexcept { return loci_; }
    inline double h2_gen() const noexcept { return h2_gen_; }
    inline double h2_vert() const noexcept { return h2_vert_; }
    inline double h2_env() const noexcept { return h2_env_; }

    inline const double& operator()(std::size_t id, ComponentType type) const {
      if (id >= n_ind_)
        throw std::runtime_error("attempting out-of-bounds access of phenotype");
      if (type == ComponentType::GENETIC)
        return ptr_gen_[id];
      else if (type == ComponentType::ENVIRONMENTAL)
        return ptr_env_[id];
      return ptr_vert_[id];
    }

    inline const double& operator()(std::size_t id) const {
      if (id >= n_ind_)
        throw std::runtime_error("attempting out-of-bounds access of phenotype");
      return vals_[id];
    }

    inline void transmit_vert() {
      double scale = std::sqrt(h2_vert_);
      std::copy(vals_.begin(), vals_.end(), ptr_vert_);
      for (std::size_t ind = 0; ind < n_ind_; ind++)
        vals_[ind] *= scale;
    }

    inline void score_values() {
      for (std::size_t ind = 0; ind < n_ind_; ind++)
        vals_[ind] = ptr_gen_[ind] + ptr_env_[ind] + ptr_vert_[ind];
    }

    void score_bitwise(Genome& genome);
    void score_tiled64(Genome& genome);
    void score(Genome& genome);

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

    std::vector<double> vals_;
  };
}
