#include <numeric>
#include <vector>
#include <cmath>

#include <amsim/stats.h>

#include <gtest/gtest.h>

double tol = 1e-12;

TEST(Stats, TestSum) {
  std::vector<double> x = {1.0, 2.0, 3.0, 4.0};

  double sum = amsim::stats::sum(x.size(), x.data(), 1);

  ASSERT_EQ(sum, 10.0);
}

TEST(Stats, TestMean) {
  std::vector<double> x = {1.0, 2.0, 3.0, 4.0};

  double mean = amsim::stats::mean(x.size(), x.data(), 1);

  ASSERT_EQ(mean, 2.5);
}

TEST(Stats, TestVar) {
  std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};

  // bessel-corrected (sample) variance
  double sample_var = amsim::stats::var(x.size(), x.data(), 1, false);

  // population variance
  double population_var = amsim::stats::var(x.size(), x.data(), 1, true);

  ASSERT_EQ(sample_var, 2.5);
  ASSERT_EQ(population_var, 2.0);
}

TEST(Stats, TestCentre) {
  std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> v_opt(x.size());
  std::vector<double> v_cnt(x.size());
  double centre_x = amsim::stats::mean(x.size(), x.data(), 1);

  amsim::stats::centre(x.size(), x.data(), 1, v_opt.data(), 1);
  amsim::stats::centre(x.size(), x.data(), 1, v_cnt.data(), 1, centre_x);

  double mean_v_opt = amsim::stats::mean(v_opt.size(), v_opt.data(), 1);
  double mean_v_centre = amsim::stats::mean(v_cnt.size(), v_cnt.data(), 1);

  ASSERT_NEAR(mean_v_opt, 0.0, tol);
  ASSERT_NEAR(mean_v_centre, 0.0, tol);
}

TEST(Stats, TestScale) {
  std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> v(x.size());

  double sd_x = std::sqrt(amsim::stats::var(x.size(), x.data(), 1));

  amsim::stats::scale(x.size(), x.data(), 1, v.data(), 1);
  amsim::stats::scale(x.size(), x.data(), 1, v.data(), 1, sd_x);

  double scale_v_opt = amsim::stats::var(v.size(), v.data(), 1);
  double scale_v_sd = amsim::stats::var(v.size(), v.data(), 1);

  ASSERT_NEAR(scale_v_opt, 1.0, tol);
  ASSERT_NEAR(scale_v_sd, 1.0, tol);
}

TEST(Stats, TestStandardise) {
  std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
  amsim::stats::standardise(x.size(), x.data(), 1, x.data(), 1);

  ASSERT_NEAR(amsim::stats::mean(x.size(), x.data(), 1), 0.0, tol);
  ASSERT_NEAR(amsim::stats::var(x.size(), x.data(), 1), 1.0, tol);
}

TEST(Stats, CorVectorIdentical) {
  std::vector<double> x = {2, 4, 6, 8};

  double r = amsim::stats::cor(x.size(), x.data(), 1, x.data(), 1);

  ASSERT_NEAR(r, 1.0, tol);
}

TEST(Stats, CorVectorWithCentreScale) {
  std::vector<double> x = {1, 2, 3, 4};
  std::vector<double> y = {2, 4, 6, 8};

  double mx = amsim::stats::mean(x.size(), x.data(), 1);
  double my = amsim::stats::mean(y.size(), y.data(), 1);
  double sx = std::sqrt(amsim::stats::var(x.size(), x.data(), 1, true));
  double sy = std::sqrt(amsim::stats::var(y.size(), y.data(), 1, true));

  double r =
      amsim::stats::cor(x.size(), x.data(), 1, y.data(), 1, mx, sx, my, sy);

  ASSERT_NEAR(r, 1.0, tol);
}

TEST(Stats, CorVectorConstantX) {
  std::vector<double> x = {1, 1, 1, 1};
  std::vector<double> y = {1, 2, 3, 4};

  double r = amsim::stats::cor(x.size(), x.data(), 1, y.data(), 1);

  ASSERT_TRUE(std::isnan(r));
}

TEST(Stats, CorMatrixCross) {
  int n = 4;
  int p = 2;
  int q = 2;

  double x[4 * 2] = {1, 2, 3, 4, 2, 3, 4, 5};

  double y[4 * 2] = {2, 4, 6, 8, 1, 2, 3, 4};

  double r[2 * 2];

  amsim::stats::cor(n, p, q, x, n, y, n, r, p);

  double r00 = r[0 + (0 * p)];
  double r01 = r[0 + (1 * p)];
  double r10 = r[1 + (0 * p)];
  double r11 = r[1 + (1 * p)];

  ASSERT_NEAR(r00, 1.0, tol);
  ASSERT_NEAR(r11, 1.0, tol);
  ASSERT_NEAR(r01, r10, tol);
  ASSERT_LE(std::abs(r01), 1.0 + tol);
}

TEST(Stats, CorMatrixSelf) {
  int n = 4;
  int p = 2;

  double x[4 * 2] = {1, 2, 3, 4, 2, 3, 4, 5};

  double r[2 * 2];

  amsim::stats::cor(n, p, x, n, r, p);

  ASSERT_NEAR(r[0 + (0 * p)], 1.0, tol);
  ASSERT_NEAR(r[1 + (1 * p)], 1.0, tol);
  ASSERT_NEAR(r[0 + (1 * p)], 1.0, tol);
  ASSERT_NEAR(r[1 + (0 * p)], 1.0, tol);
}

TEST(Stats, QuantileTest) {
  std::vector<double> x(10);
  std::iota(x.begin(), x.end(), 1);

  double q025 = amsim::stats::quantile(0.025, 10, x.data(), 1);
  double q500 = amsim::stats::quantile(0.500, 10, x.data(), 1);
  double q975 = amsim::stats::quantile(0.975, 10, x.data(), 1);

  ASSERT_NEAR(q025, 1.225, tol);
  ASSERT_NEAR(q500, 5.500, tol);
  ASSERT_NEAR(q975, 9.775, tol);
}
