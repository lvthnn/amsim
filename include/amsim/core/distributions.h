#pragma once

#include <amsim/core/rng.h>

#include <Eigen/Dense>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/math/distributions.hpp>
#include <functional>

namespace amsim {

template <typename D>
Eigen::VectorXd generate_dist(const D& dist, std::size_t n) {
  Eigen::VectorXd random(n);
  rng::UniformRange::fill(random.data(), n);
  std::ranges::transform(random, random.begin(), [dist](double u) {
    return boost::math::quantile(dist, u);
  });
  return random;
}

class RademacherDistribution {
 public:
  explicit RademacherDistribution() = default;

  static Eigen::VectorXd generate(std::size_t n) {
    Eigen::VectorXd random(n);
    rng::NormalPolar::fill(random.data(), n);
    random = random.array().sign();
    return random;
  }
};

class UniformDistribution {
 public:
  explicit UniformDistribution(double lo = 0, double hi = 1)
      : lo_(lo), hi_(hi) {}

  Eigen::VectorXd generate(std::size_t n) const {
    Eigen::VectorXd random(n);
    rng::UniformRange::fill(random.data(), n);
    random = ((hi_ - lo_) * random).array() - lo_;
    return random;
  }

 private:
  double lo_;
  double hi_;
};

// slightly redundant, but ok
class NormalDistribution {
 public:
  explicit NormalDistribution(double mean = 0, double stddev = 1)
      : mean_(mean), stddev_(stddev) {};

  Eigen::VectorXd generate(std::size_t n) const {
    Eigen::VectorXd random(n);
    rng::NormalPolar::fill(random.data(), n);
    random = mean_ + (stddev_ * random).array();
    return random;
  }

 private:
  double mean_;
  double stddev_;
};

class BetaDistribution {
 public:
  explicit BetaDistribution(double alpha = 1, double beta = 1)
      : dist_(alpha, beta) {}

  Eigen::VectorXd generate(std::size_t n) const {
    return rng::generate_dist(dist_, n);
  }

 private:
  boost::math::beta_distribution<double> dist_;
};

class ExponentialDistribution {
 public:
  explicit ExponentialDistribution(double lambda) : dist_(lambda) {}

  Eigen::VectorXd generate(std::size_t n) const {
    return generate_dist(dist_, n);
  }

 private:
  boost::math::exponential_distribution<double> dist_;
};

class GammaDistribution {
 public:
  explicit GammaDistribution(double shape, double scale)
      : dist_(shape, scale) {}

  Eigen::VectorXd generate(std::size_t n) const {
    return generate_dist(dist_, n);
  }

 private:
  boost::math::gamma_distribution<double> dist_;
};

class LaplaceDistribution {
 public:
  explicit LaplaceDistribution(double location, double scale)
      : dist_(location, scale) {}

  Eigen::VectorXd generate(std::size_t n) const {
    return generate_dist(dist_, n);
  }

 private:
  boost::math::laplace_distribution<double> dist_;
};

class StudentsTDistribution {
 public:
  explicit StudentsTDistribution(double df) : dist_(df) {}

  Eigen::VectorXd generate(std::size_t n) const {
    return generate_dist(dist_, n);
  }

 private:
  boost::math::students_t_distribution<double> dist_;
};

using Distribution = std::function<Eigen::VectorXd(std::size_t)>;

inline Distribution make_distribution(
    const std::string& dist_name,
    const std::vector<double>& params,
    bool is_probability = false) {
  // does the random variable assume values in the unit interval?
  if (is_probability) {
    std::vector<std::string> valid = {"uniform", "beta"};
    if (std::ranges::find(valid, dist_name) == valid.end())
      throw std::runtime_error(
          "Distributions for [0,1]-supported random variables must be one of "
          "'uniform(0,1)' or 'beta(a,b)'");
    if (dist_name == "uniform" && (params[0] < 0 || params[1] > 1 ||
        params[0] >= params[1]))
      throw std::runtime_error(
          "[0,1]-supported random variable with uniform distribution must have "
          "0 <= lo < hi <= 1");
  }
  if (dist_name == "rademacher") {
    if (params.size() != 0)
      throw std::runtime_error("Rademacher distribution has no parameters.");
    return [](std::size_t n) {
      return amsim::RademacherDistribution::generate(n);
    };
  }
  if (dist_name == "uniform") {
    if (params.size() != 2)
      throw std::runtime_error("Uniform distribution requires two parameters");
    return [dist = UniformDistribution(params[0], params[1])](std::size_t n) {
      return dist.generate(n);
    };
  }
  if (dist_name == "beta") {
    if (params.size() != 2)
      throw std::runtime_error("Beta distribution requires two parameters");
    return [dist = BetaDistribution(params[0], params[1])](std::size_t n) {
      return dist.generate(n);
    };
  }
  if (dist_name == "exponential") {
    if (params.size() != 1)
      throw std::runtime_error(
          "Exponential distribution requires one parameter");
    return [dist = ExponentialDistribution(params[0])](std::size_t n) {
      return dist.generate(n);
    };
  }
  if (dist_name == "gamma") {
    if (params.size() != 2)
      throw std::runtime_error("Gamma distribution requires two parameters");
    return [dist = GammaDistribution(params[0], params[1])](std::size_t n) {
      return dist.generate(n);
    };
  }
  if (dist_name == "normal") {
    if (params.size() != 2)
      throw std::runtime_error("Normal distribution requires two parameters");
    if (params[1] <= 0)
      throw std::runtime_error(
          "Standard deviance must be positive in normal distribution");
    return [dist = NormalDistribution(params[0], params[1])](std::size_t n) {
      return dist.generate(n);
    };
  }
  if (dist_name == "laplace") {
    if (params.size() != 2)
      throw std::runtime_error("Laplace distribution requires two parameters");
    return [dist = LaplaceDistribution(params[0], params[1])](std::size_t n) {
      return dist.generate(n);
    };
  }
  if (dist_name == "students_t") {
    if (params.size() != 1)
      throw std::runtime_error(
          "Student's t distribution requires one parameter");
    return [dist = StudentsTDistribution(params[0])](std::size_t n) {
      return dist.generate(n);
    };
  }
  throw std::runtime_error(std::format("Unknown distribution {}", dist_name));
}

}  // namespace amsim
