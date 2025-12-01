#include <amsim/simulation_results.h>

#include <iostream>

int main() {
  const std::filesystem::path out_dir = "amsim_multithread";
  amsim::SimulationResults results(out_dir);
  results.summarise();

  amsim::ResultsTable pheno_h2 = results("pheno_h2");
  amsim::ResultsTable pheno_gen_cor = results("pheno_gen_cor");

  results.save(
      std::vector<std::string>({"pheno_h2", "pheno_gen_cor", "pheno_tot_xcor"}),
      "amsim_results");
}
