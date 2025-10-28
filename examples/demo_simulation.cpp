#include <cstddef>
#include <vector>

#include <amsim/mating.h>
#include <amsim/simulation_builder.h>
#include <amsim/component_type.h>
#include <amsim/metricspec.h>

int main() {
  std::vector<double> v_maf(8000, 0.5);
  std::vector<double> v_rec(8000, 0.5);
  std::vector<double> v_mut(8000, 1e-8);

  amsim::SimulationBuilder builder;

  builder.simulation(10, 32000, "amsim_1610", 1234568ull)
         .genome(8000, v_maf, v_rec, v_mut)
         .phenome(2,
                  {"height", "weight"},
                  {4000, 4000},
                  {0.5, 0.5},
                  {0.5, 0.5},
                  {0.0, 0.0},
                  {1.0, 0.5,
                   0.5, 1.0},
                  {1.0, 0.5,
                   0.5, 1.0})
          .mating(amsim::MatingType::ASSORTATIVE, 2e6, 0.5, 0.9999,
                  std::vector<double>{0.4, 0.3,  // [ 0.4  0.2
                                      0.2, 0.5}) //   0.3  0.5 ]
          .metrics({
            amsim::pheno_h2(),
            amsim::pheno_comp_cor(amsim::ComponentType::GENETIC),
            amsim::pheno_comp_cor(amsim::ComponentType::ENVIRONMENTAL),
            amsim::pheno_comp_xcor(amsim::ComponentType::GENETIC),
            amsim::pheno_latent_h2(),
            amsim::pheno_latent_comp_cor(amsim::ComponentType::GENETIC),
            amsim::pheno_latent_comp_xcor(amsim::ComponentType::GENETIC)
          });

  amsim::Simulation simulation = builder.build();

  simulation.run();
}
