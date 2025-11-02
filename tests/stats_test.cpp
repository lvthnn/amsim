#include <gtest/gtest.h>

#include <vector>

#include <amsim/stats.h>

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

  // centre not passed vs passed 
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

  // scale not passed vs passed
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

  double r = amsim::stats::cor(x.size(), x.data(), 1, y.data(), 1, mx, sx, my, sy);

  ASSERT_NEAR(r, 1.0, tol);
}

TEST(Stats, CorVectorConstantX) {
  std::vector<double> x = {1, 1, 1, 1};
  std::vector<double> y = {1, 2, 3, 4};

  double r = amsim::stats::cor(x.size(), x.data(), 1, y.data(), 1);

  ASSERT_TRUE(std::isnan(r));
}

TEST(Stats, CorMatrixCross) {
  int N = 4, P = 2, Q = 2;

  double X[4*2] = {
    // col0
    1, 2, 3, 4,
    // col1
    2, 3, 4, 5
  };

  double Y[4*2] = {
    // col0 = 2 * col0 of X
    2, 4, 6, 8,
    // col1 = col1 of X - 1
    1, 2, 3, 4
  };

  double R[2*2];

  amsim::stats::cor(N, P, Q, X, N, Y, N, R, P);

  double r00 = R[0 + 0*P];
  double r01 = R[0 + 1*P];
  double r10 = R[1 + 0*P];
  double r11 = R[1 + 1*P];

  ASSERT_NEAR(r00, 1.0, tol);
  ASSERT_NEAR(r11, 1.0, tol);
  ASSERT_NEAR(r01, r10, tol);
  ASSERT_LE(std::abs(r01), 1.0 + tol);
}

TEST(Stats, CorMatrixSelf) {
  int N = 4, P = 2;

  double X[4*2] = {
    1, 2, 3, 4,
    2, 3, 4, 5
  };

  double R[2*2];

  amsim::stats::cor(N, P, X, N, R, P);

  ASSERT_NEAR(R[0 + 0*P], 1.0, tol);
  ASSERT_NEAR(R[1 + 1*P], 1.0, tol);
  ASSERT_NEAR(R[0 + 1*P], 1.0, tol);
  ASSERT_NEAR(R[1 + 0*P], 1.0, tol);
}
