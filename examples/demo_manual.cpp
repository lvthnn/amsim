#include <amsim/simulation.h>
#include <iostream>
#include <vector>
#include <filesystem>

int main() {
    std::size_t n_gen   = 2;
    std::size_t n_ind   = 2000;      // divisible by 2
    std::size_t n_loc   = 1000;
    std::size_t n_pheno = 2;

    std::vector<double> v_maf(n_loc, 0.5);
    std::vector<double> v_rec(n_loc, 0.5);
    std::vector<double> v_mut(n_loc, 1e-8);

    std::vector<std::string> v_name = {"height", "weight"};
    std::vector<std::size_t> v_n_loc = {500, 500};
    std::vector<double> v_h2_gen  = {0.5, 0.5};
    std::vector<double> v_h2_env  = {0.5, 0.5};
    std::vector<double> v_h2_vert = {0.0, 0.0};  // no VC for now

    // positive semi-definite 2x2
    std::vector<double> gen_cor = {1.0, 0.3,
                                   0.3, 1.0};
    std::vector<double> env_cor = {1.0, 0.2,
                                   0.2, 1.0};

    // assortative mating matrix 2x2
    std::vector<double> mate_cor = {0.4, 0.2,
                                    0.2, 0.5};

    std::size_t n_itr   = 1e5;      // annealing steps
    double temp_init    = 1.0;
    double temp_decay   = 0.9999;

    bool require_latent = false;

    std::vector<amsim::MetricSpec> specs = {
        amsim::pheno_h2(),
        amsim::pheno_comp_cor(amsim::ComponentType::GENETIC)
    };

    std::filesystem::path out_dir = "manual_sim_output";

    try {
        amsim::Simulation sim(
            n_gen,
            n_ind,
            out_dir,
            1234567ULL,     // seed
            n_loc,
            v_maf,
            v_rec,
            v_mut,
            n_pheno,
            v_name,
            v_n_loc,
            v_h2_gen,
            v_h2_env,
            v_h2_vert,
            gen_cor,
            env_cor,
            mate_cor,
            n_itr,
            temp_init,
            temp_decay,
            specs,
            require_latent
        );

        std::cerr << "Running manual constructor test…\n";
        sim.run();
        std::cerr << "Done.\n";

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return 1;
    }

    // Crude smoke test for output existence
    for (auto &entry : std::filesystem::directory_iterator(out_dir)) {
        std::cout << "Created: " << entry.path() << "\n";
    }
}
