#include <amsim/simulation_config.h>
#include <amsim/metricspec.h>
#include <gtest/gtest.h>

#include <vector>

class SimulationInitConfig : public ::testing::Test {
 protected:
  amsim::SimulationConfig config;
};

class SimulationGenomeConfig : public ::testing::Test {
 protected:
  void SetUp() override {
    config.simulation(100, 10000, "/tmp/test", 12345ull);
  }

  amsim::SimulationConfig config;
};

class SimulationPhenomeConfig : public ::testing::Test {
 protected:
  void SetUp() override {
    config.simulation(100, 10000, "/tmp/test", 12345ull);
    config.genome(
        100,
        std::vector<double>(100, 0.5),
        std::vector<double>(100, 0.01),
        std::vector<double>(100, 1e-5));
  }

  amsim::SimulationConfig config;
};

class SimulationMatingConfig : public ::testing::Test {
 protected:
  void SetUp() override {
    config.simulation(100, 10000, "/tmp/test", 12345ull);
    config.genome(
        100,
        std::vector<double>(100, 0.5),
        std::vector<double>(100, 0.01),
        std::vector<double>(100, 1e-5));
    config.phenome(
        2,
        {"height", "weight"},
        {50, 50},
        {0.5, 0.6},
        {0.3, 0.2},
        {0.2, 0.2},
        {1.0, 0.3, 0.3, 1.0},
        {1.0, 0.2, 0.2, 1.0});
  }

  amsim::SimulationConfig config;
};

// Test SimulationConfig.simulation

TEST_F(SimulationInitConfig, Simulation_ValidParameters) {
  EXPECT_NO_THROW(config.simulation(100, 1000, "/tmp/test", 12345ull));
  EXPECT_EQ(config.n_gen, 100);
  EXPECT_EQ(config.n_ind, 1000);
  EXPECT_EQ(config.out_dir, std::filesystem::path("/tmp/test"));
  EXPECT_EQ(config.rng_seed, 12345ull);
}

TEST_F(SimulationInitConfig, Simulation_RoundsUpOddPopulation) {
  config.simulation(1, 99, "/tmp/test", 123ull);
  EXPECT_EQ(config.n_ind, 100);
}

TEST_F(SimulationInitConfig, Simulation_Chaining) {
  auto& ref = config.simulation(1, 100, "/tmp/test", 1);
  EXPECT_EQ(&ref, &config);
}

// Test SimulationConfig.genome

TEST_F(SimulationGenomeConfig, Genome_ValidParameters) {
  EXPECT_NO_THROW(
      config.genome(3, {0.2, 0.5, 0.3}, {0.1, 0.5, 0.3}, {0.9, 0.7, 0.2}));
}

TEST_F(SimulationGenomeConfig, Genome_MismatchedVectorSizes) {
  EXPECT_THROW(
      config.genome(3, {0.2, 0.5}, {0.1, 0.5, 0.3}, {0.9, 0.7, 0.2}),
      std::invalid_argument);

  EXPECT_THROW(
      config.genome(3, {0.2, 0.5, 0.3}, {0.3}, {0.9, 0.7, 0.2}),
      std::invalid_argument);

  EXPECT_THROW(
      config.genome(3, {0.2, 0.5, 0.3}, {0.3}, {0.9, 0.7}),
      std::invalid_argument);
}

TEST_F(SimulationGenomeConfig, Genome_InvalidMAF) {
  EXPECT_THROW(config.genome(1, {-0.5}, {0.5}, {0.5}), std::invalid_argument);
  EXPECT_THROW(config.genome(1, {1.5}, {0.5}, {0.5}), std::invalid_argument);
}

TEST_F(SimulationGenomeConfig, Genome_InvalidRecombination) {
  EXPECT_THROW(config.genome(1, {0.5}, {-0.5}, {0.5}), std::invalid_argument);
  EXPECT_THROW(config.genome(1, {0.5}, {1.5}, {0.5}), std::invalid_argument);
}

TEST_F(SimulationGenomeConfig, Genome_InvalidMutation) {
  EXPECT_THROW(config.genome(1, {0.5}, {0.5}, {-0.5}), std::invalid_argument);
  EXPECT_THROW(config.genome(1, {0.5}, {0.5}, {1.5}), std::invalid_argument);
}

TEST_F(SimulationGenomeConfig, Genome_Chaining) {
  auto& ref = config.genome(1, {0.5}, {0.5}, {0.5});
  EXPECT_EQ(&ref, &config);
}

// Test SimulationConfig.phenome

TEST_F(SimulationPhenomeConfig, Phenome_ValidParameters) {
  EXPECT_NO_THROW(config.phenome(
      2,
      {"height", "weight"},
      {50, 50},
      {0.5, 0.6},
      {0.3, 0.2},
      {0.2, 0.2},
      {1.0, 0.3, 0.3, 1.0},
      {1.0, 0.2, 0.2, 1.0}));
}

TEST_F(SimulationPhenomeConfig, Phenome_WrongNumberOfNames) {
  EXPECT_THROW(
      config.phenome(
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

TEST_F(SimulationPhenomeConfig, Phenome_WrongNumberOfLoci) {
  EXPECT_THROW(
      config.phenome(
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

TEST_F(SimulationPhenomeConfig, Phenome_WrongNumberOfH2Gen) {
  EXPECT_THROW(
      config.phenome(
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

TEST_F(SimulationPhenomeConfig, Phenome_WrongNumberOfH2Env) {
  EXPECT_THROW(
      config.phenome(
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

TEST_F(SimulationPhenomeConfig, Phenome_WrongNumberOfH2Vert) {
  EXPECT_THROW(
      config.phenome(
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

TEST_F(SimulationPhenomeConfig, Phenome_TooManyCausalLoci) {
  EXPECT_THROW(
      config.phenome(1, {"height"}, {200}, {0.5}, {0.3}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, Phenome_InvalidH2Gen) {
  EXPECT_THROW(
      config.phenome(1, {"height"}, {50}, {1.5}, {0.3}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
  EXPECT_THROW(
      config.phenome(1, {"height"}, {50}, {-0.5}, {0.3}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, Phenome_InvalidH2Env) {
  EXPECT_THROW(
      config.phenome(1, {"height"}, {50}, {0.5}, {1.5}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
  EXPECT_THROW(
      config.phenome(1, {"height"}, {50}, {0.5}, {-0.5}, {0.2}, {1.0}, {1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, Phenome_InvalidH2Vert) {
  EXPECT_THROW(
      config.phenome(1, {"height"}, {50}, {0.5}, {0.3}, {1.5}, {1.0}, {1.0}),
      std::invalid_argument);
  EXPECT_THROW(
      config.phenome(1, {"height"}, {50}, {0.5}, {0.3}, {-0.5}, {1.0}, {1.0}),
      std::invalid_argument);
}

TEST_F(SimulationPhenomeConfig, Phenome_InvalidGenCor) {
  EXPECT_THROW(
      config.phenome(
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

TEST_F(SimulationPhenomeConfig, Phenome_InvalidEnvCor) {
  EXPECT_THROW(
      config.phenome(
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

TEST_F(SimulationPhenomeConfig, Phenome_Chaining) {
  auto& ref =
      config.phenome(1, {"height"}, {50}, {0.5}, {0.3}, {0.2}, {1.0}, {1.0});
  EXPECT_EQ(&ref, &config);
}

// Test SimulationConfig.mating

TEST_F(SimulationMatingConfig, Mating_RandomMating) {
  EXPECT_NO_THROW(config.mating(amsim::MatingType::RANDOM, std::nullopt,
                                 std::nullopt, std::nullopt, std::nullopt));
  EXPECT_EQ(config.n_itr, 0);
  EXPECT_EQ(config.mate_cor.size(), 4);
}

TEST_F(SimulationMatingConfig, Mating_AssortativeWithoutCorrelation) {
  EXPECT_THROW(config.mating(amsim::MatingType::ASSORTATIVE, 100, 1.0, 0.95,
                             std::nullopt),
               std::runtime_error);
}

TEST_F(SimulationMatingConfig, Mating_AssortativeWithValidCorrelation) {
  EXPECT_NO_THROW(config.mating(amsim::MatingType::ASSORTATIVE, 100, 1.0, 0.95,
                                 std::vector<double>{0.2, 0.5, 0.3, 0.4}));
}

TEST_F(SimulationMatingConfig, Mating_NegativeTemperature) {
  EXPECT_THROW(config.mating(amsim::MatingType::ASSORTATIVE, 100, -1.0, 0.95,
                             std::vector<double>{0.2, 0.5, 0.3, 0.4}),
               std::invalid_argument);
}

TEST_F(SimulationMatingConfig, Mating_NegativeDecay) {
  EXPECT_THROW(config.mating(amsim::MatingType::ASSORTATIVE, 100, 1.0, -0.95,
                             std::vector<double>{0.2, 0.5, 0.3, 0.4}),
               std::invalid_argument);
}

TEST_F(SimulationMatingConfig, Mating_InfeasibleCorrelation) {
  EXPECT_THROW(config.mating(amsim::MatingType::ASSORTATIVE, 100, 1.0, 0.95,
                             std::vector<double>{0.99, 0.99, 0.99, 0.99}),
               std::invalid_argument);
}

TEST_F(SimulationMatingConfig, Mating_Chaining) {
  auto& ref = config.mating(amsim::MatingType::RANDOM, std::nullopt,
                            std::nullopt, std::nullopt, std::nullopt);
  EXPECT_EQ(&ref, &config);
}

// Test SimulationConfig.metrics

TEST_F(SimulationMatingConfig, Metrics_EmptyList) {
  config.metrics({});
  EXPECT_FALSE(config.require_lat);
}

TEST_F(SimulationMatingConfig, Metrics_WithLatentRequirement) {
  amsim::MetricSpec spec = amsim::pheno_latent_h2();
  config.metrics({spec});
  EXPECT_TRUE(config.require_lat);
}

TEST_F(SimulationMatingConfig, Metrics_Chaining) {
  auto& ref = config.metrics({});
  EXPECT_EQ(&ref, &config);
}
