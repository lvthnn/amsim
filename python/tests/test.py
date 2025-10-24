from pathlib import Path

from amsim import Simulation, MatingType, ComponentType
from amsim import metrics 

if __name__ == '__main__':
    result_dir = Path('/Users/karihlynsson/Documents/amsimpy_test')

    # configure the simulation
    simulation = Simulation().simulation(
        n_generations=10,
        n_individuals=256000,
        output_dir=result_dir,
        random_seed=12345671284124888
    ).genome(
        n_loci=20000,
        locus_mafs=0.5,
        locus_recombination=0.5,
        locus_mutation=1e-8
    ).phenome(
        n_phenotypes=5,
        names=['y1', 'y2', 'y3', 'y4', 'y5'],
        loci=4000,
        h2_genetic=0.5,
        h2_environmental=0.5,
        h2_vertical=0.0,
        # 5x5
        genetic_cor= [
            [1.0, 0.8, 0.6, 0.4, 0.2],
            [0.8, 1.0, 0.5, 0.3, 0.1],
            [0.6, 0.5, 1.0, 0.2, 0.0],
            [0.4, 0.3, 0.2, 1.0, 0.4],
            [0.2, 0.1, 0.0, 0.4, 1.0]
        ],
        environmental_cor= [
            [1.0, 0.2, 0.1, 0.0, 0.0],
            [0.2, 1.0, 0.3, 0.1, 0.0],
            [0.1, 0.3, 1.0, 0.2, 0.1],
            [0.0, 0.1, 0.2, 1.0, 0.3],
            [0.0, 0.0, 0.1, 0.3, 1.0]
        ]
    ).mating(
        mating_type=MatingType.ASSORTATIVE,
        n_iterations=2_000_000,
        temp_init=0.5,
        temp_decay=0.9999,
        mate_cor = [ # cross-correlation matrix, no unit diagonal
            [0.0, 0.3, 0.2, 0.1, 0.0],
            [0.3, 0.0, 0.3, 0.2, 0.1],
            [0.2, 0.3, 0.0, 0.3, 0.2],
            [0.1, 0.2, 0.3, 0.0, 0.3],
            [0.0, 0.1, 0.2, 0.3, 0.0]
        ]
    ).metrics(
        metrics=[
            metrics.pheno_h2(),
            metrics.pheno_comp_cor(ComponentType.GENETIC),
            metrics.pheno_comp_cor(ComponentType.ENVIRONMENTAL),
            metrics.pheno_comp_xcor(ComponentType.GENETIC),
            metrics.pheno_comp_xcor(ComponentType.TOTAL),
            metrics.pheno_latent_h2(),
            metrics.pheno_latent_comp_cor(ComponentType.GENETIC),
            metrics.pheno_latent_comp_cor(ComponentType.ENVIRONMENTAL),
            metrics.pheno_latent_comp_xcor(ComponentType.GENETIC),
            metrics.pheno_latent_comp_xcor(ComponentType.TOTAL)
        ]
    )

    simulation.run(n_replicates=10, n_threads=10)
