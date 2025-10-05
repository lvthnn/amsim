#pragma once
#include <vector>
#include <cstddef>
#include <random>
#include <tuple>

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

  class GeneralModel : public MatingModel {
  public:
    GeneralModel(std::vector<std::vector<double> const*> vals_ptr,
                 std::tuple<std::size_t, std::size_t, double> cor_,
                 const std::size_t n_itr, const std::size_t n_sex,
                 double tmp_init, double tmp_decay);

    void update_vals(std::vector<std::vector<double> const*> vals_ptr);
    std::vector<std::size_t> match() override;

  private:
    MatingType type_ = MatingType::ASSORTATIVE;
    std::vector<std::vector<double> const*> vals_ptr_;
    std::tuple<std::size_t, std::size_t, double> cor_;
    std::mt19937 g_;
    const std::size_t n_itr_;
    const std::size_t n_sex_;
    double tmp_init_;
    double tmp_decay_;

    double delta_(std::vector<std::size_t> cur, std::size_t i0, std::size_t i1);
    double energy_(std::vector<std::size_t> matching);
  };
}
