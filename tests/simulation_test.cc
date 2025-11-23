#include <amsim/metricspec.h>
#include <amsim/simulation_config.h>
#include <gtest/gtest.h>

#include <vector>

class SimulationInitConfig : public ::testing::Test {
 protected:
  amsim::SimulationConfig config_;
};

class SimulationGenomeConfig : public ::testing::Test {
 protected:
  void SetUp() override {
    config_.simulation(100, 10000, "/tmp/test", 12345ULL);
  }

  amsim::SimulationConfig config_;
};

class SimulationPhenomeConfig : public ::testing::Test {
 protected:
  void SetUp() override {
    config_.simulation(100, 10000, "/tmp/test", 12345ULL);
    config_.genome(
        100,
        std::vector<double>(100, 0.5),
        std::vector<double>(100, 0.01),
        std::vector<double>(100, 1e-5));
  }

  amsim::SimulationConfig config_;
};

class SimulationMatingConfig : public ::testing::Test {
 protected:
  void SetUp() override {
    config_.simulation(100, 10000, "/tmp/test", 12345ULL);
    config_.genome(
        100,
        std::vector<double>(100, 0.5),
        std::vector<double>(100, 0.01),
        std::vector<double>(100, 1e-5));
    config_.phenome(
        2,
        {"height", "weight"},
        {50, 50},
        {0.5, 0.6},
        {0.3, 0.2},
        {0.2, 0.2},
        {1.0, 0.3, 0.3, 1.0},
        {1.0, 0.2, 0.2, 1.0});
  }

  amsim::SimulationConfig config_;
};

// Test SimulationConfig.simulation

TEST_F(SimulationInitConfig, SimulationValidParametrs) {
  EXPECT_NO_THROW(config_.simulation(100, 1000, "/tmp/test", 12345ULL));
  EXPECT_EQ(config_.n_gen, 100);
  EXPECT_EQ(config_.n_ind, 1000);
  EXPECT_EQ(config_.out_dir, std::filesystem::path("/tmp/test"));
  EXPECT_EQ(config_.rng_seed, 12345ULL);
}

TEST_F(SimulationInitConfig, SimulationRoundsUpOddPopulation) {
  config_.simulation(1, 99, "/tmp/test", 123ULL);
  EXPECT_EQ(config_.n_ind, 100);
}

TEST_F(SimulationInitConfig, SimulationChaining) {
  auto& ref = config_.simulation(1, 100, "/tmp/test", 1);
  EXPECT_EQ(&ref, &config_);
}

// Test SimulationConfig.genome

TEST_F(SimulationGenomeConfig, GenomeValidParameters) {
  EXPECT_NO_THROW(
      config_.genome(3, {0.2, 0.5, 0.3}, {0.1, 0.5, 0.3}, {0.9, 0.7, 0.2}));
}

TEST_F(SimulationGenomeConfig, GenomeMismatchedVectorSizes) {
  EXPECT_THROW(
      config_.genome(3, {0.2, 0.5}, {0.1, 0.5, 0.3}, {0.9, 0.7, 0.2}),
      std::invalid_argument);

  EXPECT_THROW(
      config_.genome(3, {0.2, 0.5, 0.3}, {0.3}, {0.9, 0.7, 0.2}),
      std::invalid_argument);

  EXPECT_THROW(
      config_.genome(3, {0.2, 0.5, 0.3}, {0.3}, {0.9, 0.7}),
      std::invalid_argument);
}

TEST_F(SimulationGenomeConfig, GenomeInvalidMAF) {
  EXPECT_THROW(config_.genome(1, {-0.5}, {0.5}, {0.5}), std::invalid_argument);
  EXPECT_THROW(config_.genome(1, {1.5}, {0.5}, {0.5}), std::invalid_argument);
}

TEST_F(SimulationGenomeConfig, GenomeInvalidRecombination) {
  EXPECT_THROW(config_.genome(1, {0.5}, {-0.5}, {0.5}), std::invalid_argument);
  EXPECT_THROW(config_.genome(1, {0.5}, {1.5}, {0.5}), std::invalid_argument);
}

TEST_F(SimulationGenomeConfig, GenomeInvalidMutation) {
  EXPECT_THROW(config_.genome(1, {0.5}, {0.5}, {-0.5}), std::invalid_argument);
  EXPECT_THROW(config_.genome(1, {0.5}, {0.5}, {1.5}), std::invalid_argument);
}

TEST_F(SimulationGenomeConfig, GenomeChaining) {
  auto& ref = config_.genome(1, {0.5}, {0.5}, {0.5});
  EXPECT_EQ(&ref, &config_);
}

// Test SimulationConfig.phenome

TEST_F(SimulationPhenomeConfig, PhenomeValidParameters) {
  EXPECT_NO_THROW(config_.phenome(
      2,
      {"height", "weight"},
      {50, 50},
      {0.5, 0.6},
      {0.3, 0.2},
      {0.2, 0.2},
      {1.0, 0.3, 0.3, 1.0},
      {1.0, 0.2, 0.2, 1.0}));
}

TEST_F(SimulationPhenomeConfig, PhenomeWrongNumberOfNames) {
  EXPECT_THROW(
      config_.phenome(
          2,
          {"height"},
          {50, 50},
          {0.5, 0.6},
          {0.3, 0.2},
          {0.2, 0.2},
          {1.0, 0.3, 0.3, 1.0},
          {1.0, 0.2, 0.2, 1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeWrongNumberOfLoci) {
  EXPECT_THROW(
      config_.phenome(
          2,
          {"height", "weight"},
          {50},
          {0.5, 0.6},
          {0.3, 0.2},
          {0.2, 0.2},
          {1.0, 0.3, 0.3, 1.0},
          {1.0, 0.2, 0.2, 1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeWrongNumberOfH2Gen) {
  EXPECT_THROW(
      config_.phenome(
          2,
          {"height", "weight"},
          {50, 50},
          {0.5},
          {0.3, 0.2},
          {0.2, 0.2},
          {1.0, 0.3, 0.3, 1.0},
          {1.0, 0.2, 0.2, 1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeWrongNumberOfH2Env) {
  EXPECT_THROW(
      config_.phenome(
          2,
          {"height", "weight"},
          {50, 50},
          {0.5, 0.6},
          {0.3},
          {0.2, 0.2},
          {1.0, 0.3, 0.3, 1.0},
          {1.0, 0.2, 0.2, 1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeWrongNumberOfH2Vert) {
  EXPECT_THROW(
      config_.phenome(
          2,
          {"height", "weight"},
          {50, 50},
          {0.5, 0.6},
          {0.3, 0.2},
          {0.2},
          {1.0, 0.3, 0.3, 1.0},
          {1.0, 0.2, 0.2, 1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeTooManyCausalLoci) {
  EXPECT_THROW(
      config_.phenome(1, {"height"}, {200}, {0.5}, {0.3}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeInvalidH2Gen) {
  EXPECT_THROW(
      config_.phenome(1, {"height"}, {50}, {1.5}, {0.3}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
  EXPECT_THROW(
      config_.phenome(1, {"height"}, {50}, {-0.5}, {0.3}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeInvalidH2Env) {
  EXPECT_THROW(
      config_.phenome(1, {"height"}, {50}, {0.5}, {1.5}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
  EXPECT_THROW(
      config_.phenome(1, {"height"}, {50}, {0.5}, {-0.5}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeInvalidH2Vert) {
  EXPECT_THROW(
      config_.phenome(1, {"height"}, {50}, {0.5}, {0.3}, {1.5}, {1.0}, {1.0}),
      std::invalid_argument);
  EXPECT_THROW(
      config_.phenome(1, {"height"}, {50}, {0.5}, {0.3}, {-0.5}, {1.0}, {1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeInvalidGenCor) {
  EXPECT_THROW(
      config_.phenome(
          2,
          {"height", "weight"},
          {50, 50},
          {0.5, 0.6},
          {0.3, 0.2},
          {0.2, 0.2},
          {1.0, 1.5, 1.5, 1.0},
          {1.0, 0.2, 0.2, 1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeInvalidEnvCor) {
  EXPECT_THROW(
      config_.phenome(
          2,
          {"height", "weight"},
          {50, 50},
          {0.5, 0.6},
          {0.3, 0.2},
          {0.2, 0.2},
          {1.0, 0.3, 0.3, 1.0},
          {1.0, -1.5, -1.5, 1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, PhenomeChaining) {
  auto& ref =
      config_.phenome(1, {"height"}, {50}, {0.5}, {0.3}, {0.2}, {1.0}, {1.0});
  EXPECT_EQ(&ref, &config_);
}

// Test SimulationConfig.mating

TEST_F(SimulationMatingConfig, MatingRandomMating) {
  EXPECT_NO_THROW(config_.random_mating());
  EXPECT_EQ(config_.n_itr, 0);
  EXPECT_EQ(config_.mate_cor.size(), 4);
}

TEST_F(SimulationMatingConfig, MatingAssortativeWithValidCorrelation) {
  EXPECT_NO_THROW(config_.assortative_mating(
      std::vector<double>{0.2, 0.5, 0.3, 0.4},
      1e-6,
      100,
      1.0,
      0.95));
}

TEST_F(SimulationMatingConfig, MatingInfeasibleCorrelation) {
  EXPECT_THROW(
      config_.assortative_mating(
          std::vector<double>{0.99, 0.99, 0.99, 0.99},
          1e-6,
          100,
          1.0,
          0.95),
      std::invalid_argument);
}

TEST_F(SimulationMatingConfig, MatingNegativeTemperature) {
  EXPECT_THROW(
      config_.assortative_mating(
          std::vector<double>{0.2, 0.5, 0.3, 0.4},
          1e-6,
          100,
          -1.0,
          0.95),
      std::invalid_argument);
}

TEST_F(SimulationMatingConfig, MatingNegativeDecay) {
  EXPECT_THROW(
      config_.assortative_mating(
          std::vector<double>{0.2, 0.5, 0.3, 0.4},
          1e-6,
          100,
          1.0,
          -0.95),
      std::invalid_argument);
}

TEST_F(SimulationMatingConfig, MatingChaining) {
  auto& ref = config_.random_mating();
  EXPECT_EQ(&ref, &config_);
}

// Test SimulationConfig.metrics

TEST_F(SimulationMatingConfig, MetricsEmptyList) {
  config_.metrics({});
  EXPECT_FALSE(config_.require_lat);
}

TEST_F(SimulationMatingConfig, MetricsWithLatentRequirement) {
  amsim::MetricSpec spec = amsim::pheno_latent_h2();
  config_.metrics({spec});
  EXPECT_TRUE(config_.require_lat);
}

TEST_F(SimulationMatingConfig, MetricsChaining) {
  auto& ref = config_.metrics({});
  EXPECT_EQ(&ref, &config_);
}
