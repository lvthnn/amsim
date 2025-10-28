#ifndef AMSIMCPP_SIMULATION_H
#define AMSIMCPP_SIMULATION_H

#pragma once
#include <vector>
#include <fstream>
#include <filesystem>

#include <amsim/genome.h>
#include <amsim/phenobuf.h>
#include <amsim/phenoarch.h>
#include <amsim/phenobuf.h>
#include <amsim/mating.h>
#include <amsim/metricspec.h>
#include <amsim/mating.h>
#include <amsim/simulation_status.h>

namespace amsim {

  class Simulation {
  public:
    Simulation(// SIMULATION PARAMETERS
							 std::size_t              n_gen,
							 std::size_t              n_ind,
							 std::filesystem::path    out_dir,
							 std::uint64_t            rng_seed,
							 // GENOME PARAMETERS
							 std::size_t              n_loc,
							 std::vector<double>      v_maf,
							 std::vector<double>      v_rec,
							 std::vector<double>      v_mut,
							 // PHENOME PARAMETERS
							 std::size_t              n_pheno,
							 std::vector<std::string> v_name,
							 std::vector<std::size_t> v_n_loc,
							 std::vector<double>      v_h2_gen,
							 std::vector<double>      v_h2_env,
						   std::vector<double>      v_h2_vert,
							 std::vector<double>      gen_cor,
							 std::vector<double>      env_cor,
							 // MATING MODEL PARAMETERS
               std::vector<double>      mate_cor,
							 std::size_t              n_itr,
							 double                   tmp_init,
							 double                   tmp_decay,
							 std::vector<MetricSpec>  specs,
							 bool                     require_latent);

    Simulation(Simulation&& other) noexcept
      : n_gen(other.n_gen),
        out_dir(std::move(other.out_dir)),
        status_(other.status_),
        rng_(std::move(other.rng_)),
        genome_(std::move(other.genome_)),
        arch_(std::move(other.arch_)),
        buf_(std::move(other.buf_)),
        phenotypes_(std::move(other.phenotypes_)),
        model_(std::move(other.model_)),
        ctx_(genome_, arch_, buf_, phenotypes_, model_),
        metrics_(std::move(other.metrics_)),
        streams_(std::move(other.streams_))
    {}

    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation& operator=(Simulation&&) = delete;

    const std::size_t n_gen;
    const std::filesystem::path out_dir;

    void run();

  private:
    SimulationStatus    status_;
    rng::Xoshiro256ss   rng_;
    Genome              genome_;
    PhenoArch           arch_;
    PhenoBuf            buf_;
    PhenotypeList       phenotypes_;
    AssortativeModel    model_;
    SimulationContext   ctx_;
    std::vector<Metric> metrics_;

    std::vector<std::unique_ptr<std::ofstream>> streams_;
    void stream_(std::size_t gen);
  };

}

#endif //AMSIMCPP_SIMULATION_H
