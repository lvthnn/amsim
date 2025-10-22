from pathlib import Path

from amsim import SimulationBuilder, MatingType, ComponentType
from amsim import metrics 

if __name__ == '__main__':
    result_dir = Path('/Users/modulor/Documents/amsimpy_test')

    # initialise the simulation builder
    builder = SimulationBuilder()

    # configure the simulation
    builder.simulation(
        n_generations=10,
        n_individuals=32000,
        output_dir=result_dir,
        random_seed=12345671284124888
    ).genome(
        n_loci=8000,
        locus_mafs=0.5,
        locus_recombination=0.5,
        locus_mutation=1e-8
    ).phenome(
        n_phenotypes=2,
        names=["height", "weight"],
        loci=4000,
        h2_genetic=0.5,
        h2_environmental=0.5,
        h2_vertical=0.0,
        genetic_cor=[[1.0, 0.5], [0.5, 1.0]],
        environmental_cor=[[1.0, 0.5], [0.5, 1.0]]
    ).mating(
        mating_type=MatingType.ASSORTATIVE,
        n_iterations=2_000_000,
        temp_init=0.5,
        temp_decay=0.9999,
        mate_cor=[[0.4, 0.3], [0.2, 0.5]]
    ).metrics(
        metric_specs=[
            metrics.pheno_h2(2),
            metrics.pheno_comp_cor(2, ComponentType.GENETIC),
            metrics.pheno_comp_cor(2, ComponentType.ENVIRONMENTAL),
            metrics.pheno_comp_xcor(2, ComponentType.TOTAL),
            metrics.pheno_latent_h2(2),
            metrics.pheno_latent_comp_cor(2, ComponentType.GENETIC),
            metrics.pheno_latent_comp_xcor(2, ComponentType.GENETIC)
        ]
    )

    # build the simulation and run it
    simulation = builder.build().run(n_replicates=10, n_threads=5)
