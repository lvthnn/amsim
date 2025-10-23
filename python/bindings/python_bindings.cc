#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <amsim/simulation.h>
#include <amsim/simulation_builder.h>
#include <amsim/component_type.h>
#include <amsim/metric.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
  // enums required for setup
  py::enum_<amsim::MatingType>(m, "MatingType")
    .value("RANDOM", amsim::MatingType::RANDOM)
    .value("ASSORTATIVE", amsim::MatingType::ASSORTATIVE)
    .export_values();

  py::enum_<amsim::ComponentType>(m, "ComponentType")
    .value("GENETIC", amsim::ComponentType::GENETIC)
    .value("ENVIRONMENTAL", amsim::ComponentType::ENVIRONMENTAL)
    .value("VERTICAL", amsim::ComponentType::VERTICAL)
    .value("TOTAL", amsim::ComponentType::TOTAL)
    .export_values();

  // metrics used by builder
  py::class_<amsim::Metric>(m, "_Metric")
    .def_readonly("name", &amsim::Metric::name);

  m.def("_loc_maf", &amsim::metrics::loc_maf,
        py::arg("n_loci"),
        "Create locus MAF metric");

  m.def("_loc_mean", &amsim::metrics::loc_mean,
        py::arg("n_loci"),
        "Create locus mean metric");

  m.def("_loc_var", &amsim::metrics::loc_var,
        py::arg("n_loci"),
        "Create locus variance metric");

  m.def("_pheno_h2", &amsim::metrics::pheno_h2,
        py::arg("n_phenotypes"),
        "Create heritability metric"); 

  m.def("_pheno_comp_mean", &amsim::metrics::pheno_comp_mean,
        py::arg("n_phenotypes"),
        py::arg("component_type"),
        "Create component mean metric");

  m.def("_pheno_comp_var", &amsim::metrics::pheno_comp_var,
        py::arg("n_phenotypes"),
        py::arg("component_type"),
        "Create component variance metric");

  m.def("_pheno_comp_cor", &amsim::metrics::pheno_comp_cor,
        py::arg("n_phenotypes"),
        py::arg("component_type"),
        "Create component correlation matrix metric");

  m.def("_pheno_comp_xcor", &amsim::metrics::pheno_comp_xcor,
        py::arg("n_phenotypes"),
        py::arg("component_type"),
        "Create between-mate component correlation matrix metric");

  m.def("_pheno_latent_h2", &amsim::metrics::pheno_latent_h2,
        py::arg("n_phenotypes"),
        "Create latent phenotype heritability metric");

  m.def("_pheno_latent_comp_cor", &amsim::metrics::pheno_latent_comp_cor,
        py::arg("n_phenotypes"),
        py::arg("component_type"),
        "Create latent component correlation matrix metric");

  m.def("_pheno_latent_comp_xcor", &amsim::metrics::pheno_latent_comp_xcor,
        py::arg("n_phenotypes"),
        py::arg("component_type"),
        "Create between-mate latent component correlation matrix metric");

  // simulation builder and simulation classes
  py::class_<amsim::SimulationBuilder>(m, "_SimulationBuilder")
    .def(py::init<>())
    .def("simulation", &amsim::SimulationBuilder::simulation,
         py::arg("n_generations"),
         py::arg("n_individuals"),
         py::arg("output_dir"),
         py::arg("random_seed"),
         py::return_value_policy::reference_internal,
         R"doc(
             Configure simulation parameters
            
             Parameters
             ----------
             n_generations : int
                Number of generations to simulate
             n_individuals : int
                Number of individuals per generation
             output_dir : str
                Directory for output files
             random_seed : int
                Random seed for reproducibility

             Returns
             -------
             SimulationBuilder
                 Self, for method chaining
         )doc")
    .def("genome", &amsim::SimulationBuilder::genome,
         py::arg("n_loci"),
         py::arg("locus_mafs"),
         py::arg("locus_recombination"),
         py::arg("locus_mutation"),
         py::return_value_policy::reference_internal,
         R"doc(
             Configure genome parameters
            
             Parameters
             ----------
             n_loci : int
                Number of genetic loci to simulate
             locus_mafs : list of float
                Locus minor allele frequency vector 
             locus_recombination : list of float 
                Locus recombination probabilities 
             locus_mutation : list of float
                Locus mutation probabilities

             Returns
             -------
             SimulationBuilder
                 Self, for method chaining
         )doc")
    .def("phenome", &amsim::SimulationBuilder::phenome,
         py::arg("n_phenotypes"),
         py::arg("names"),
         py::arg("loci"),
         py::arg("h2_genetic"),
         py::arg("h2_environmental"),
         py::arg("h2_vertical"),
         py::arg("genetic_cor"),
         py::arg("environmental_cor"),
         py::return_value_policy::reference_internal)
    .def("mating", &amsim::SimulationBuilder::mating,
         py::arg("mating_type"),
         py::arg("n_iterations"),
         py::arg("temp_init"),
         py::arg("temp_decay"),
         py::arg("mate_cor"),
         py::return_value_policy::reference_internal)
    .def("metrics", &amsim::SimulationBuilder::metrics,
         py::arg("metrics"),
         py::return_value_policy::reference_internal)
    .def("build", &amsim::SimulationBuilder::build);

  py::class_<amsim::Simulation>(m, "_Simulation")
    .def("run", &amsim::Simulation::run);
}
