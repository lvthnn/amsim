#pragma once

#include <amsim/state.h>
#include <amsim/params.h>

#include <filesystem>

namespace amsim {

void setup_output(std::filesystem::path out_dir);

// Preprocess parameters for initialisers
void simulation_preprocess(Params& params);

// Perform a single simulation transformation step
void simulation_run(
    State& state,
    const Params& params,
    std::size_t n_gen,
    std::filesystem::path out_dir);

void simulation_single();

void simulation_multithread();

}
