#include <iostream>

#include <amsim/simulation_results.h>

int main() {
  const std::filesystem::path out_dir = "amsim_multithread";
  amsim::SimulationResults results(out_dir);

  // std::cout << "n_replicates: " << results.n_replicates_ << "\n";
  // for (const auto& rep_dir : results.rep_dirs_) {
  //   std::cout << rep_dir.string() << "\n";
  // }

  // for (const auto& metric_name : results.metric_names_) {
  //   std::cout << metric_name << "\n";
  // }

  // std::cout << results.metric_names_.size();

  // results.summarise_metric_("pheno_h2");
  results.summarise();
  results.print("pheno_tot_xcor");
}
