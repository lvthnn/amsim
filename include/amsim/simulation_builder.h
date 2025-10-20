#ifndef AMSIMCPP_SIMULATION_BUILDER_H
#define AMSIMCPP_SIMULATION_BUILDER_H

#include <amsim/rng.h>
#include <amsim/simulation.h>
#include <amsim/simulation_status.h>
#include <amsim/haplobuf.h>
#include <amsim/genome.h>
#include <amsim/phenotype.h>
#include <amsim/phenobuf.h>
#include <amsim/phenoarch.h>
#include <amsim/mating.h>
#include <amsim/metric.h>

namespace amsim {

  class SimulationBuilder {
  public:
    SimulationBuilder();

    SimulationBuilder& simulation(std::size_t n_gen, std::size_t n_ind,
                                  std::string out_dir, std::uint64_t rng_seed);

    SimulationBuilder& genome(std::size_t n_loc, std::vector<double> v_maf,
                              std::vector<double> v_rec,
                              std::vector<double> v_mut);

    SimulationBuilder& phenome(std::size_t n_pheno,
                               std::vector<std::string> v_name,
                               std::vector<std::size_t> v_n_loc,
                               std::vector<double> v_h2_gen,
                               std::vector<double> v_h2_env,
                               std::vector<double> v_h2_vert,
                               std::vector<double> gen_cor,
                               std::vector<double> env_cor);

    SimulationBuilder& mating(MatingType type, std::optional<std::size_t> n_itr,
                              std::optional<double> tmp_init,
                              std::optional<double> tmp_decay,
                              std::optional<std::vector<double>> mate_cor);

    SimulationBuilder& metrics(std::vector<Metric> metrics);

    Simulation build();

  private:
    //-- TRACK SIMULATION PROGRESS ------------------------------------------//
    SimulationStatus         status_;      // simulation progress
    //-----------------------------------------------------------------------//

    //-- SIMULATION PARAMETERS ----------------------------------------------//
    std::size_t              n_gen_;       // number of generations simulated
    std::size_t              n_ind_;       // number of individuals
    std::filesystem::path    out_dir_;     // name of folder to write data to
    std::uint64_t            rng_seed_;    // seed of xoshiro random device
    //-----------------------------------------------------------------------//

    //-- GENOME PARAMETERS --------------------------------------------------//
    std::size_t              n_loc_;       // number of loci to simulate
    std::vector<double>      v_maf_;       // locus mafs
    std::vector<double>      v_rec_;       // recombination map
    std::vector<double>      v_mut_;       // mutation map
    //-----------------------------------------------------------------------//

    //-- PHENOTYPE PARAMETERS -----------------------------------------------//
    std::size_t              n_pheno_;     // number of phenotypes
    std::vector<std::string> v_name_;      // phenotype names
    std::vector<std::size_t> v_n_loc_;     // number of loci per phenotype
    std::vector<double>      v_h2_gen_;    // genetic component variance
    std::vector<double>      v_h2_env_;    // environmental component variance
    std::vector<double>      v_h2_vert_;   // vertical component variance
    std::vector<double>      gen_cor_;     // genetic component corr.
    std::vector<double>      env_cor_;     // environmental component corr.
    //-----------------------------------------------------------------------//

    //-- MATING MODEL PARAMETERS --------------------------------------------//
    std::size_t              n_itr_;       // number of optimisation steps
    double                   tmp_init_;    // initial annealing temperature
    double                   tmp_decay_;   // annealing temperature decay
    std::vector<double>      mate_cor_;    // mate correlation matrix
    //-----------------------------------------------------------------------//

    //-- METRIC PARMETERS ---------------------------------------------------//
    std::vector<Metric>      metrics_;
    bool                     require_lat_ = false;
    //-----------------------------------------------------------------------//
  };

}

#endif // AMSIMCPP_SIMULATION_BUILDER_H
