#ifndef AMSIMCPP_MATING_H
#define AMSIMCPP_MATING_H

#pragma once
#include <vector>
#include <cstddef>
#include <random>

#include <amsim/rng.h>
#include <amsim/phenotype.h>

namespace amsim {
  enum class MatingType {
    RANDOM,
    ASSORTATIVE
  };

  class MatingModel {
  public:
    explicit MatingModel(MatingType type, std::size_t n_sex)
      : type_(type), g_(std::random_device{}()), n_sex_(n_sex) {}

    virtual ~MatingModel() = default;
    virtual std::vector<std::size_t> match() = 0;

  protected:
    MatingType type_;
    std::mt19937 g_;
    std::size_t n_sex_;

    std::vector<std::size_t> rand_state_();
  };

  class AssortativeModel : public MatingModel {
  public:
    AssortativeModel(const PhenotypeList &phenotypes,
                     std::vector<double> cor, const std::size_t n_itr,
                     const std::size_t n_sex, const rng::Xoshiro256ss &rng,
                     double tmp_init = 1e-9, double tmp_decay = 0.99999999);

    void display_cor();
    std::vector<std::size_t> match() override;
    void update(const PhenotypeList &phenotypes);

  private:
    std::vector<const double*> ptr_tot_;

    const std::vector<double> cor_;
    const std::size_t n_pheno_;
    const std::size_t n_sex_;
    const std::size_t n_itr_;
    const double tmp_init_;
    const double tmp_decay_;

    std::vector<double> male_;
    std::vector<double> female_;
    std::vector<std::size_t> state_;

    void arrange_();

    std::vector<double> compute_cor_();
    std::vector<double> compute_delta_(std::size_t i0, std::size_t i1);

    double compute_denergy_(const std::vector<double> &cur,
                            const std::vector<double> &target,
                            const std::vector<double> &delta);

    rng::UniformIntRange swap_;
    rng::UniformRange acc_;
  };
}
#endif
