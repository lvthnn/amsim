#pragma once

#include <amsim/core/params.h>

#include <Eigen/Dense>
#include <cstddef>
#include <vector>

namespace amsim {

enum Component { Genetic = 0, Environmental = 1, Nurture = 2, Total = 3 };

inline Component operator++(Component& type, int) {
  Component old = type;
  type = (type == Component::Total)
             ? Component::Genetic
             : Component(static_cast<int>(type) + 1);
  return old;
}

inline Component& operator++(Component& type) {
  type = (type == Component::Total)
             ? Component::Genetic
             : Component(static_cast<int>(type) + 1);
  return type;
}

inline std::string to_string(Component type) {
  switch (type) {
    case Component::Genetic:
      return "genetic";
    case Component::Environmental:
      return "environ";
    case Component::Nurture:
      return "nurture";
    case Component::Total:
      return "total";
  }
  __builtin_unreachable();
}

inline std::ostream& operator<<(std::ostream& os, Component type) {
  return os << to_string(type);
}

// Buffer structure to store phenotype data
class PhenoBuf {
 public:
  explicit PhenoBuf(const Params& params);

  Eigen::Map<const Eigen::MatrixXd> operator()() const {
    const double* pos = data_.data();
    return Eigen::Map<const Eigen::MatrixXd>(pos, n_ind_, n_pheno_ * 4);
  }

  Eigen::Map<Eigen::MatrixXd> operator()() {
    double* pos = data_.data();
    return Eigen::Map<Eigen::MatrixXd>(pos, n_ind_, n_pheno_ * 4);
  }

  // Retrieve const matrix of values for specified component
  Eigen::Map<const Eigen::MatrixXd> operator()(Component type = Component::Total) const {
    const double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type))];
    return Eigen::Map<const Eigen::MatrixXd>(pos, n_ind_, n_pheno_);
  }

  // Retrieve matrix of values for specified component
  Eigen::Map<Eigen::MatrixXd> operator()(Component type = Component::Total) {
    double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type))];
    return Eigen::Map<Eigen::MatrixXd>(pos, n_ind_, n_pheno_);
  }

  // Retrieve const vector of values for specified phenotype component
  Eigen::Map<const Eigen::VectorXd> operator()(
      std::size_t id, Component type = Component::Total) const {
    const double* pos =
        &data_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    return Eigen::Map<const Eigen::VectorXd>(pos, n_ind_);
  }

  // Retrieve vector of values for specified phenotype component
  Eigen::Map<Eigen::VectorXd> operator()(std::size_t id, Component type = Component::Total) {
    double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    return Eigen::Map<Eigen::VectorXd>(pos, n_ind_);
  }

  auto male(Component type = Component::Total) {
    Eigen::OuterStride<> stride(n_ind_);
    double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type))];
    return Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(
        pos, n_sex_, n_pheno_, stride);
  }

  auto female(Component type = Component::Total) {
    Eigen::OuterStride<> stride(n_ind_);
    double* pos =
        &data_[(n_ind_ * (n_pheno_ * static_cast<int>(type))) + (n_ind_ / 2)];
    return Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(
        pos, n_sex_, n_pheno_, stride);
  }

  Eigen::Map<Eigen::VectorXd> male(
      std::size_t id, Component type = Component::Total) {
    double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    return Eigen::Map<Eigen::VectorXd>(pos, n_sex_);
  }

  Eigen::Map<Eigen::VectorXd> female(
      std::size_t id, Component type = Component::Total) {
    double* pos =
        &data_[(n_ind_ * (n_pheno_ * static_cast<int>(type) + id)) + n_sex_];
    return Eigen::Map<Eigen::VectorXd>(pos, n_sex_);
  }

  auto male(Component type = Component::Total) const {
    Eigen::OuterStride<> stride(n_ind_);
    const double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type))];
    return Eigen::Map<const Eigen::MatrixXd, 0, Eigen::OuterStride<>>(
        pos, n_sex_, n_pheno_, stride);
  }

  auto female(Component type = Component::Total) const {
    Eigen::OuterStride<> stride(n_ind_);
    const double* pos =
        &data_[(n_ind_ * (n_pheno_ * static_cast<int>(type))) + (n_ind_ / 2)];
    return Eigen::Map<const Eigen::MatrixXd, 0, Eigen::OuterStride<>>(
        pos, n_sex_, n_pheno_, stride);
  }

  Eigen::Map<const Eigen::VectorXd> male(
      std::size_t id, Component type = Component::Total) const {
    const double* pos =
        &data_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    return Eigen::Map<const Eigen::VectorXd>(pos, n_sex_);
  }

  Eigen::Map<const Eigen::VectorXd> female(
      std::size_t id, Component type = Component::Total) const {
    const double* pos =
        &data_[(n_ind_ * (n_pheno_ * static_cast<int>(type) + id)) + n_sex_];
    return Eigen::Map<const Eigen::VectorXd>(pos, n_sex_);
  }

  // means and variances of phenotype components
  void compute_stats();

  double comp_mean(std::size_t pheno_id, Component type) const {
    return comp_mean_(static_cast<int>(type), pheno_id);
  }

  double comp_var(std::size_t pheno_id, Component type) const {
    return comp_var_(static_cast<int>(type), pheno_id);
  }

  Eigen::Ref<const Eigen::VectorXd> comp_mean(Component type) const {
    return comp_mean_.row(static_cast<int>(type));
  }

  Eigen::Ref<const Eigen::VectorXd> comp_var(Component type) const {
    return comp_var_.row(static_cast<int>(type));
  }

 private:
  const std::size_t n_ind_;  ///< Number of individuals
  const std::size_t n_sex_;
  const std::size_t n_pheno_;  ///< Number of phenotypes
  std::vector<double> data_;   ///< Main phenotype component buffer

  Eigen::MatrixXd comp_mean_;
  Eigen::MatrixXd comp_var_;
};

}  // namespace amsim
