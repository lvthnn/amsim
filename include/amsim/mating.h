#ifndef AMSIMCPP_MATING_H
#define AMSIMCPP_MATING_H

#pragma once
#include <vector>
#include <cstddef>
#include <random>

#include <amsim/phenotype.h>

namespace amsim::mating {
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
    AssortativeModel(const std::vector<double*> &ptr_tot, std::vector<double> cor,
                 const std::size_t n_itr, const std::size_t n_sex,
                 double tmp_init = 1e-9, double tmp_decay = 0.99995);

    void update_vals(std::vector<std::vector<double> const*> ptr_tot);
    std::vector<std::size_t> match() override;

  private:
    const std::vector<double> cor_;
    const std::vector<double*> ptr_tot_;
    const std::size_t n_pheno_;
    const std::size_t n_sex_;
    const std::size_t n_itr_;
    const double tmp_init_;
    const double tmp_decay_;

    std::vector<double> male_;
    std::vector<double> female_;
    std::vector<std::size_t> opt_state_;

    std::vector<double> cmp_cor_();
    void setup_(std::vector<std::size_t> state);
    double delta_(std::vector<std::size_t> cur, std::size_t i0, std::size_t i1);
  };
}
#endif