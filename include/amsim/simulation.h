#pragma once

#include <amsim/output/estimator.h>
#include <amsim/params.h>
#include <amsim/state.h>

#include <filesystem>

namespace amsim {

void setup_output(std::filesystem::path out_dir);

// Preprocess parameters for initialisers
void simulation_preprocess(Params& params);

// Perform a single simulation transformation step
void simulation_run(
    State& state,
    const Params& params,
    const Estimators& estimators,
    std::size_t n_gen,
    std::optional<std::size_t> rep_id = std::nullopt);

void simulation_single();

void simulation_multithread();

}  // namespace amsim
