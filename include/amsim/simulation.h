#pragma once

#include <amsim/output/estimator.h>
#include <amsim/params.h>
#include <amsim/state.h>
#include <amsim/setup.h>

#include <filesystem>

namespace amsim {

void setup_output(std::filesystem::path out_dir);

// Build and run a simulation
void simulation_build(Simulation& simulation);

// Preprocess parameters for initialisers
void simulation_preprocess(Params& params);

// Perform a single simulation transformation step
void simulation_run(
    State& state,
    Params& params,
    const Estimators& estimators,
    std::size_t n_gen,
    std::optional<std::size_t> rep_id = std::nullopt);

void simulation_single();

void simulation_multithread();

}  // namespace amsim
