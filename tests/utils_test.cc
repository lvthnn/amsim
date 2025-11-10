#include <amsim/utils.h>
#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

// Tests for amsim::utils::assert_probs

TEST(Utils, AssertProbs_ValidProbabilities) {
  std::vector<double> probs = {0.0, 0.5, 1.0, 0.25, 0.75};
  EXPECT_NO_THROW(amsim::utils::assert_probs(probs.size(), probs.data(), 1));
}

TEST(Utils, AssertProbs_AllZeros) {
  std::vector<double> probs = {0.0, 0.0, 0.0};
  EXPECT_NO_THROW(amsim::utils::assert_probs(probs.size(), probs.data(), 1));
}

TEST(Utils, AssertProbs_AllOnes) {
  std::vector<double> probs = {1.0, 1.0, 1.0};
  EXPECT_NO_THROW(amsim::utils::assert_probs(probs.size(), probs.data(), 1));
}

TEST(Utils, AssertProbs_EmptyArray) {
  std::vector<double> probs;
  EXPECT_NO_THROW(amsim::utils::assert_probs(0, probs.data(), 1));
}

TEST(Utils, AssertProbs_SingleElement) {
  double prob = 0.5;
  EXPECT_NO_THROW(amsim::utils::assert_probs(1, &prob, 1));
}

TEST(Utils, AssertProbs_NegativeValue) {
  std::vector<double> probs = {0.5, -0.1, 0.3};
  EXPECT_THROW(
      {
        try {
          amsim::utils::assert_probs(probs.size(), probs.data(), 1);
        } catch (const std::invalid_argument& e) {
          EXPECT_STREQ("invalid element 1; must be in [0, 1]", e.what());
          throw;
        }
      },
      std::invalid_argument);
}

TEST(Utils, AssertProbs_ValueGreaterThanOne) {
  std::vector<double> probs = {0.2, 0.8, 1.1};
  EXPECT_THROW(
      {
        try {
          amsim::utils::assert_probs(probs.size(), probs.data(), 1);
        } catch (const std::invalid_argument& e) {
          EXPECT_STREQ("invalid element 2; must be in [0, 1]", e.what());
          throw;
        }
      },
      std::invalid_argument);
}

TEST(Utils, AssertProbs_FirstElementInvalid) {
  std::vector<double> probs = {-0.5, 0.5, 0.8};
  EXPECT_THROW(
      {
        try {
          amsim::utils::assert_probs(probs.size(), probs.data(), 1);
        } catch (const std::invalid_argument& e) {
          EXPECT_STREQ("invalid element 0; must be in [0, 1]", e.what());
          throw;
        }
      },
      std::invalid_argument);
}

TEST(Utils, AssertProbs_WithStride) {
  std::vector<double> data = {0.3, 999.0, 0.5, 999.0, 0.8};
  EXPECT_NO_THROW(amsim::utils::assert_probs(3, data.data(), 2));
}

TEST(Utils, AssertProbs_WithStrideInvalid) {
  std::vector<double> data = {0.3, 999.0, 1.5, 999.0, 0.8};
  EXPECT_THROW(
      amsim::utils::assert_probs(3, data.data(), 2), std::invalid_argument);
}

// Tests for amsim::utils::assert_cors

TEST(Utils, AssertCors_ValidCorrelations) {
  std::vector<double> cors = {-1.0, -0.5, 0.0, 0.5, 1.0, 0.75, -0.25};
  EXPECT_NO_THROW(amsim::utils::assert_cors(cors.size(), cors.data(), 1));
}

TEST(Utils, AssertCors_AllZeros) {
  std::vector<double> cors = {0.0, 0.0, 0.0};
  EXPECT_NO_THROW(amsim::utils::assert_cors(cors.size(), cors.data(), 1));
}

TEST(Utils, AssertCors_AllOnes) {
  std::vector<double> cors = {1.0, 1.0, 1.0};
  EXPECT_NO_THROW(amsim::utils::assert_cors(cors.size(), cors.data(), 1));
}

TEST(Utils, AssertCors_AllNegativeOnes) {
  std::vector<double> cors = {-1.0, -1.0, -1.0};
  EXPECT_NO_THROW(amsim::utils::assert_cors(cors.size(), cors.data(), 1));
}

TEST(Utils, AssertCors_EmptyArray) {
  std::vector<double> cors;
  EXPECT_NO_THROW(amsim::utils::assert_cors(0, cors.data(), 1));
}

TEST(Utils, AssertCors_SingleElement) {
  double cor = 0.5;
  EXPECT_NO_THROW(amsim::utils::assert_cors(1, &cor, 1));
}

TEST(Utils, AssertCors_ValueLessThanNegativeOne) {
  std::vector<double> cors = {0.5, -1.1, 0.3};
  EXPECT_THROW(
      {
        try {
          amsim::utils::assert_cors(cors.size(), cors.data(), 1);
        } catch (const std::invalid_argument& e) {
          EXPECT_STREQ("invalid element 1; must be in [-1, 1]", e.what());
          throw;
        }
      },
      std::invalid_argument);
}

TEST(Utils, AssertCors_ValueGreaterThanOne) {
  std::vector<double> cors = {0.2, 0.8, 1.1};
  EXPECT_THROW(
      {
        try {
          amsim::utils::assert_cors(cors.size(), cors.data(), 1);
        } catch (const std::invalid_argument& e) {
          EXPECT_STREQ("invalid element 2; must be in [-1, 1]", e.what());
          throw;
        }
      },
      std::invalid_argument);
}

TEST(Utils, AssertCors_FirstElementInvalid) {
  std::vector<double> cors = {-1.5, 0.5, 0.8};
  EXPECT_THROW(
      {
        try {
          amsim::utils::assert_cors(cors.size(), cors.data(), 1);
        } catch (const std::invalid_argument& e) {
          EXPECT_STREQ("invalid element 0; must be in [-1, 1]", e.what());
          throw;
        }
      },
      std::invalid_argument);
}

TEST(Utils, AssertCors_WithStride) {
  std::vector<double> data = {-0.8, 999.0, 0.5, 999.0, 0.9};
  EXPECT_NO_THROW(amsim::utils::assert_cors(3, data.data(), 2));
}

TEST(Utils, AssertCors_WithStrideInvalid) {
  std::vector<double> data = {0.3, 999.0, -1.5, 999.0, 0.8};
  EXPECT_THROW(
      amsim::utils::assert_cors(3, data.data(), 2), std::invalid_argument);
}

// Tests for amsim::utils::assert_udiag

TEST(Utils, AssertUDiag_IdentityMatrix) {
  std::vector<double> matrix = {1.0, 0.5, 0.3, 0.5, 1.0, 0.2, 0.3, 0.2, 1.0};
  EXPECT_NO_THROW(amsim::utils::assert_udiag(3, matrix.data(), 3));
}

TEST(Utils, AssertUDiag_SingleElement) {
  double matrix = 1.0;
  EXPECT_NO_THROW(amsim::utils::assert_udiag(1, &matrix, 1));
}

TEST(Utils, AssertUDiag_DiagonalNotOne) {
  std::vector<double> matrix = {1.0, 0.5, 0.3, 0.5, 0.9, 0.2, 0.3, 0.2, 1.0};
  EXPECT_THROW(
      {
        try {
          amsim::utils::assert_udiag(3, matrix.data(), 3);
        } catch (const std::invalid_argument& e) {
          EXPECT_STREQ("diagonal element 1 does not equal one", e.what());
          throw;
        }
      },
      std::invalid_argument);
}

TEST(Utils, AssertUDiag_FirstDiagonalWrong) {
  std::vector<double> matrix = {0.99, 0.5, 0.3, 0.5, 1.0, 0.2, 0.3, 0.2, 1.0};
  EXPECT_THROW(
      {
        try {
          amsim::utils::assert_udiag(3, matrix.data(), 3);
        } catch (const std::invalid_argument& e) {
          EXPECT_STREQ("diagonal element 0 does not equal one", e.what());
          throw;
        }
      },
      std::invalid_argument);
}

TEST(Utils, AssertUDiag_LastDiagonalWrong) {
  std::vector<double> matrix = {1.0, 0.5, 0.3, 0.5, 1.0, 0.2, 0.3, 0.2, 1.01};
  EXPECT_THROW(
      {
        try {
          amsim::utils::assert_udiag(3, matrix.data(), 3);
        } catch (const std::invalid_argument& e) {
          EXPECT_STREQ("diagonal element 2 does not equal one", e.what());
          throw;
        }
      },
      std::invalid_argument);
}

TEST(Utils, AssertUDiag_WithinTolerance) {
  std::vector<double> matrix = {
      1.0 + 1e-13, 0.5, 0.3, 0.5, 1.0 - 1e-13, 0.2, 0.3, 0.2, 1.0};
  EXPECT_NO_THROW(amsim::utils::assert_udiag(3, matrix.data(), 3));
}

TEST(Utils, AssertUDiag_OutsideTolerance) {
  std::vector<double> matrix = {
      1.0 + 2e-12, 0.5, 0.3, 0.5, 1.0, 0.2, 0.3, 0.2, 1.0};
  EXPECT_THROW(
      amsim::utils::assert_udiag(3, matrix.data(), 3), std::invalid_argument);
}

TEST(Utils, AssertUDiag_WithLeadingDimension) {
  std::vector<double> matrix = {
      1.0, 0.5, 0.3, 999.0, 0.5, 1.0, 0.2, 999.0, 0.3, 0.2, 1.0, 999.0};
  EXPECT_NO_THROW(amsim::utils::assert_udiag(3, matrix.data(), 4));
}

// Tests for amsim::utils::assert_psd

TEST(Utils, AssertPSD_IdentityMatrix) {
  std::vector<double> matrix = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  EXPECT_NO_THROW(amsim::utils::assert_psd(3, matrix.data(), 3));
}

TEST(Utils, AssertPSD_PositiveDefiniteMatrix) {
  // Symmetric positive definite matrix
  std::vector<double> matrix = {4.0, 2.0, 1.0, 2.0, 3.0, 1.0, 1.0, 1.0, 2.0};
  EXPECT_NO_THROW(amsim::utils::assert_psd(3, matrix.data(), 3));
}

TEST(Utils, AssertPSD_CorrelationMatrix) {
  // Valid correlation matrix
  std::vector<double> matrix = {1.0, 0.8, 0.5, 0.8, 1.0, 0.6, 0.5, 0.6, 1.0};
  EXPECT_NO_THROW(amsim::utils::assert_psd(3, matrix.data(), 3));
}

TEST(Utils, AssertPSD_SingleElement) {
  double matrix = 1.0;
  EXPECT_NO_THROW(amsim::utils::assert_psd(1, &matrix, 1));
}

TEST(Utils, AssertPSD_TwoByTwo) {
  std::vector<double> matrix = {1.0, 0.5, 0.5, 1.0};
  EXPECT_NO_THROW(amsim::utils::assert_psd(2, matrix.data(), 2));
}

TEST(Utils, AssertPSD_ZeroMatrix) {
  // Zero matrix is positive semi-definite (all eigenvalues are 0)
  std::vector<double> matrix(9, 0.0);
  EXPECT_NO_THROW(amsim::utils::assert_psd(3, matrix.data(), 3));
}

TEST(Utils, AssertPSD_NegativeEigenvalue) {
  // Matrix with a negative eigenvalue
  std::vector<double> matrix = {1.0, 2.0, 3.0, 2.0, -1.0, 1.0, 3.0, 1.0, -2.0};
  EXPECT_THROW(
      amsim::utils::assert_psd(3, matrix.data(), 3), std::invalid_argument);
}

TEST(Utils, AssertCrossCor_Feasible) {
  std::vector<double> matrix(16, 0.2);
  EXPECT_NO_THROW(amsim::utils::assert_cross_cor(4, matrix.data(), 4));
}

TEST(Utils, AssertCrossCor_Infeasible) {
  std::vector<double> matrix(25, 0.25);
  EXPECT_THROW(
    amsim::utils::assert_cross_cor(5, matrix.data(), 5), std::invalid_argument);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
