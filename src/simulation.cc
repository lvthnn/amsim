#include <amsim/state.h>
#include <amsim/params.h>

#include <amsim/preprocess.h>
#include <amsim/transform.h>
#include <amsim/initialise.h>

#include <memory>
#include <filesystem>

namespace amsim {

void setup_output(const std::filesystem::path& out_dir) {
  if (std::filesystem::exists(out_dir))
    throw std::runtime_error("directory " + out_dir.string() + " already exists!");
  std::filesystem::create_directory(out_dir);
}

void simulation_preprocess(Params& params) {
  preprocess::calibrate_ld_matrix(params);
  preprocess::optimise_phenotype_arch(params);
}

void simulation_run(
    State& state,
    const Params& params,
    std::size_t n_gen,
    const std::filesystem::path& out_dir) {
  // Set up the output directory
  setup_output(out_dir);

  // Founder haplotype initialiser
  std::unique_ptr<genome::HaplotypeGenerator> haplo;
  if (params.geno.has_ld()) {
    haplo = std::make_unique<genome::HaploGeneratorLD>(params); // then grapple with this demon
  } else {
    haplo = std::make_unique<genome::HaplotypeGeneratorIID>(params); // implement this first
  }

  // Set up required transformers
  phenome::PhenotypeScorer score(params); // needs to be implemented
  mating::AssortativeMating mate(params); // ready
  genome::GenomeUpdater update(params); // needs to be implemented

  // Generate the initial state
  haplo->generate_haplotypes(state.geno);

  for (std::size_t gen = 0; gen < n_gen; ++gen) {
    // Compute genome mean and variances
    state.geno.compute_mafs();
    state.geno.compute_stats();

    // Score phenotypes using operator()
    score(state.pheno);
    state.pheno.compute_stats();

    // Match mates
    auto match = mate(state.pheno.male(), state.pheno.female());

    // Update the genome
    update(state);

    // Produce subviews / sample the population
  }
}

}

