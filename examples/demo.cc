#include <cstddef>
#include <vector>
#include <chrono>
#include <random>

#include <amsim/genome.h>
#include <amsim/phenotype.h>
#include <amsim/metric.h>
#include <amsim/utils.h>
#include <amsim/mating.h>

int main() {
  // Generate random MAFs and set up seeded Xoshiro device
  std::random_device rd;
  std::mt19937 g(rd());
  std::uniform_real_distribution<double> maf_dist(0.0, 1.0);

  std::size_t rng_seed = 123442255188518518ull;
  amsim::rng::Xoshiro256ss rng = amsim::rng::seed_xoshiro(amsim::rng::auto_seed(rng_seed));

  // Parameters for genome
  std::size_t n_ind = 256000;
  std::size_t n_loc = 5000;
  std::size_t n_pheno = 3;
  std::vector<double> v_mut(n_loc, 1e-8);
  std::vector<double> v_rec(n_loc, 0.5);
  std::vector<double> v_maf(n_loc, 0.5);

  // Parameters for phenotypes
  double h2_gen  = 0.5;
  double h2_env  = 0.5;
  double h2_vert = 0.0;
  std::vector<std::size_t> v_n_loc(n_pheno, 1000);
  std::vector<double> v_h2_gen(n_pheno, h2_gen);
  std::vector<double> v_h2_env(n_pheno, h2_env);

  // Genetic component correlation (@ panmixis)
  std::vector<double> gen_cor {
    1.0, 0.2, 0.3,   // col 1
    0.2, 1.0, 0.4,   // col 2
    0.3, 0.4, 1.0    // col 3
  };

  // Environmental component correlation
  std::vector<double> env_cor {
    1.0, 0.2, 0.1,   // col 1
    0.2, 1.0, 0.1,   // col 2
    0.1, 0.1, 1.0    // col 3
  };

  // Between-mate trait correlation
  std::vector<double> mate_cor {
    0.2, 0.2, 0.2,
    0.2, 0.2, 0.2,
    0.2, 0.2, 0.2
  };

  // Run the simulation
  const auto total_start = std::chrono::steady_clock::now();

  // Declare genome, phenotype buffer and architecture
  amsim::Genome    genome(n_ind, n_loc, v_mut, v_rec, v_maf, rng_seed);
  amsim::PhenoBuf  buf(n_ind, 3);
  amsim::PhenoArch arch(n_pheno, n_loc, v_n_loc, v_h2_gen, v_h2_env, gen_cor, env_cor, rng);

  arch.optim_arch(1e4, 1e-12);

  // amsim::utils::time_step("Optimising genetic architecture", [&] { arch.optim_arch(3e3, 1e-12); });

  amsim::Phenotype height(buf, arch, "height", h2_gen, h2_env, h2_vert);
  amsim::Phenotype weight(buf, arch, "weight", h2_gen, h2_env, h2_vert);
  amsim::Phenotype bmi(buf, arch, "bmi", h2_gen, h2_env, h2_vert);

  amsim::PhenotypeList phenotypes = {height, weight, bmi};

  // initialise
  genome.generate_haplotypes();

  // typical simulation cycle: --------
  genome.compute_mafs();
  genome.compute_stats();
  arch.gen_env(buf(0, amsim::ENVIRONMENTAL), n_ind);
  height.score(genome);
  weight.score(genome);
  bmi.score(genome);
  // generate the environmental components using arch
  height.compute_stats();
  weight.compute_stats();
  bmi.compute_stats();

  // perform mating

  amsim::AssortativeModel model(phenotypes, mate_cor, 10 * n_ind, n_ind / 2, rng, 1e-9, 0.999999);
  amsim::utils::time_step("Mate matching routine", [&](){ model.match(); });
  model.display_cor();

  // stream metrics


  // update

  // -----------------------------------

  // issue: phenotype vector copies phenotypes, copying uninitialised pointers -> unexpected behaviour
  // fix this by implementing some wrapper class (?) or using std::reference_wrapper

  // for (amsim::ComponentType type = amsim::ComponentType::GENETIC; type <= amsim::ComponentType::TOTAL; type++) {
  //   std::cout << "----------------------------------------\n";
  //   std::cout << "component type: " << type << "\n";
  //   std::cout << "height (mean):  " << height.comp_mean(type) << "\n";
  //   std::cout << "weight (mean):  " << weight.comp_mean(type) << "\n";
  //   std::cout << "bmi    (mean):  " << bmi.comp_mean(type) << "\n";
  //   std::cout << "height (var):   " << height.comp_var(type) << "\n";
  //   std::cout << "weight (var):   " << weight.comp_var(type) << "\n";
  //   std::cout << "bmi    (var):   " << bmi.comp_var(type) << "\n";
  //   if (type == amsim::ComponentType::TOTAL) break;
  // }

  // amsim::Metric metric_gen_cor = amsim::metrics::comp_cor("gen_cor", n_pheno, amsim::ComponentType::GENETIC);
  // metric_gen_cor.stream(phenotypes, genome);

  // amsim::utils::time_step("Generated haplotypes", [&]{ genome.generate_haplotypes(); });
  // amsim::utils::time_step("Computed MAFs",        [&]{ genome.compute_mafs(); });
  // amsim::utils::time_step("Computed locus stats", [&]{ genome.compute_stats(); });
  // amsim::utils::time_step("Scored phenotypes",    [&]{ height.score(genome); weight.score(genome); });
  // amsim::utils::time_step("Computed phenotype stats",    [&]{ height.cmp_stats(); weight.cmp_stats(); });
  // amsim::utils::time_step("Transpose",            [&]{ genome.transpose(); });

  // const long long total_ms =
  //     std::chrono::duration_cast<std::chrono::milliseconds>(
  //         std::chrono::steady_clock::now() - total_start).count();
  // std::cout << "Total elapsed: " << total_ms << " ms\n\n";
}
