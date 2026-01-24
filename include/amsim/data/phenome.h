#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <vector>

namespace amsim::phenome {

enum ComponentType { GENETIC = 0, ENVIRONMENTAL = 1, VERTICAL = 2, TOTAL = 3 };

inline ComponentType operator++(ComponentType& type, int) {
  ComponentType old = type;
  type = (type == ComponentType::TOTAL)
             ? ComponentType::GENETIC
             : ComponentType(static_cast<int>(type) + 1);
  return old;
}

inline ComponentType& operator++(ComponentType& type) {
  type = (type == ComponentType::TOTAL)
             ? ComponentType::GENETIC
             : ComponentType(static_cast<int>(type) + 1);
  return type;
}

inline std::string to_string(ComponentType type) {
  switch (type) {
    case ComponentType::GENETIC:
      return "genetic";
    case ComponentType::ENVIRONMENTAL:
      return "environ";
    case ComponentType::VERTICAL:
      return "nurture";
    case ComponentType::TOTAL:
      return "total";
  }
  __builtin_unreachable();
}

inline std::ostream& operator<<(std::ostream& os, ComponentType type) {
  return os << to_string(type);
}

// Buffer structure to store phenotype data
class PhenoBuf {
 public:
  PhenoBuf(std::size_t n_ind, std::size_t n_pheno);

  // Retrieve const matrix of values for specified component
  Eigen::Map<const Eigen::MatrixXd> operator()(ComponentType type) const {
    const double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type))];
    return Eigen::Map<const Eigen::MatrixXd>(pos, n_ind_, n_pheno_);
  }

  // Retrieve matrix of values for specified component
  Eigen::Map<Eigen::MatrixXd> operator()(ComponentType type) {
    double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type))];
    return Eigen::Map<Eigen::MatrixXd>(pos, n_ind_, n_pheno_);
  }

  // Retrieve const vector of values for specified phenotype component
  Eigen::Map<const Eigen::VectorXd> operator()(
      std::size_t id, ComponentType type) const {
    const double* pos =
        &data_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    return Eigen::Map<const Eigen::VectorXd>(pos, n_ind_);
  }

  // Retrieve vector of values for specified phenotype component
  Eigen::Map<Eigen::VectorXd> operator()(std::size_t id, ComponentType type) {
    double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    return Eigen::Map<Eigen::VectorXd>(pos, n_ind_);
  }

  auto male(ComponentType type = ComponentType::TOTAL) {
    Eigen::OuterStride<> stride(n_ind_);
    double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type))];
    return Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(
        pos, n_sex_, n_pheno_, stride);
  }

  auto female(ComponentType type = ComponentType::TOTAL) {
    Eigen::OuterStride<> stride(n_ind_);
    double* pos =
        &data_[(n_ind_ * (n_pheno_ * static_cast<int>(type))) + (n_ind_ / 2)];
    return Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(
        pos, n_sex_, n_pheno_, stride);
  }

  Eigen::Map<Eigen::VectorXd> male(
      std::size_t id, ComponentType type = ComponentType::TOTAL) {
    double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    return Eigen::Map<Eigen::VectorXd>(pos, n_sex_);
  }

  Eigen::Map<Eigen::VectorXd> female(
      std::size_t id, ComponentType type = ComponentType::TOTAL) {
    double* pos =
        &data_[(n_ind_ * (n_pheno_ * static_cast<int>(type) + id)) + n_sex_];
    return Eigen::Map<Eigen::VectorXd>(pos, n_sex_);
  }

  auto male(ComponentType type = ComponentType::TOTAL) const {
    Eigen::OuterStride<> stride(n_ind_);
    const double* pos = &data_[n_ind_ * (n_pheno_ * static_cast<int>(type))];
    return Eigen::Map<const Eigen::MatrixXd, 0, Eigen::OuterStride<>>(
        pos, n_sex_, n_pheno_, stride);
  }

  auto female(ComponentType type = ComponentType::TOTAL) const {
    Eigen::OuterStride<> stride(n_ind_);
    const double* pos =
        &data_[(n_ind_ * (n_pheno_ * static_cast<int>(type))) + (n_ind_ / 2)];
    return Eigen::Map<const Eigen::MatrixXd, 0, Eigen::OuterStride<>>(
        pos, n_sex_, n_pheno_, stride);
  }

  Eigen::Map<const Eigen::VectorXd> male(
      std::size_t id, ComponentType type = ComponentType::TOTAL) const {
    const double* pos =
        &data_[n_ind_ * (n_pheno_ * static_cast<int>(type) + id)];
    return Eigen::Map<const Eigen::VectorXd>(pos, n_sex_);
  }

  Eigen::Map<const Eigen::VectorXd> female(
      std::size_t id, ComponentType type = ComponentType::TOTAL) const {
    const double* pos =
        &data_[(n_ind_ * (n_pheno_ * static_cast<int>(type) + id)) + n_sex_];
    return Eigen::Map<const Eigen::VectorXd>(pos, n_sex_);
  }

  // Compute means and variances of phenotype components
  void compute_stats();

 private:
  const std::size_t n_ind_;  ///< Number of individuals
  const std::size_t n_sex_;
  const std::size_t n_pheno_;  ///< Number of phenotypes
  std::vector<double> data_;   ///< Main phenotype component buffer

  Eigen::MatrixXd comp_mean_;
  Eigen::MatrixXd comp_var_;
};

}  // namespace amsim::phenome
