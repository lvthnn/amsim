#include <cstddef>
#include <vector>

#include <amsim/simulation_builder.h>

int main() {
  std::vector<double> v_maf(5000, 0.5);
  std::vector<double> v_rec(5000, 0.5);
  std::vector<double> v_mut(5000, 1e-8);

  amsim::Simulation simulation = amsim::SimulationBuilder()
    .simulation(15, 256000, "amsim_1610", 123456ull)
    .genome(5000, v_maf, v_rec, v_mut)
    .phenome(3,
             {"height", "weight", "bmi"},
             {1000, 1000, 1000},
             {0.5, 0.5, 0.5},
             {0.5, 0.5, 0.5},
             {0.0, 0.0, 0.0},
             {1.0, 0.2, 0.1,
              0.2, 1.0, 0.1,
              0.1, 0.1, 1.0},
             {1.0, 0.2, 0.1,
              0.2, 1.0, 0.1,
              0.1, 0.1, 1.0})
    .mating(amsim::MatingType::ASSORTATIVE, 1e6, 1e-9, 0.9999999,
            std::vector<double>{0.2, 0.2, 0.2,
                                0.2, 0.2, 0.2,
                                0.2, 0.2, 0.2})
    .metrics(std::vector<amsim::Metric>{ amsim::metrics::pheno_h2("pheno_heritability", 3) })
    .build();
}
