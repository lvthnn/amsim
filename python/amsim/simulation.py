"""Interface between Python and amsimcpp backend."""

from time import time_ns
from pathlib import Path
from multiprocessing import Pool
from typing import Self

import numpy as np

from amsim._core import _SimulationBuilder, _Simulation, MatingType
from amsim.utils import broadcast_values, column_major
from amsim.metrics import Metric, build_metrics


def simulate_replication(params: dict, rep_id: int | None = None) -> None:
    if rep_id is not None:
        output_dir = Path(params['simulation']['output_dir'])
        rep_dir = output_dir / f"rep_{rep_id:03d}"
        rep_dir.mkdir(parents=True, exist_ok=True)
        rep_dir  = str(rep_dir.resolve())
        rep_seed = params['simulation']['random_seed'] + rep_id
    else:
        rep_dir = params['simulation']['output_dir']
        rep_seed = params['simulation']['random_seed']

    b = _SimulationBuilder()

    b.simulation(
        n_generations=params['simulation']['n_generations'],
        n_individuals=params['simulation']['n_individuals'],
        output_dir=rep_dir,
        random_seed=rep_seed,
    )

    b.genome(
        n_loci=params['genome']['n_loci'],
        locus_mafs=broadcast_values(
            params['genome']['locus_mafs'],
            params['genome']['n_loci']
        ),
        locus_recombination=broadcast_values(
            params['genome']['locus_recombination'],
            params['genome']['n_loci']
        ),
        locus_mutation=broadcast_values(
            params['genome']['locus_mutation'],
            params['genome']['n_loci']
        ),
    )

    b.phenome(
        n_phenotypes=params['phenome']['n_phenotypes'],
        names=params['phenome']['names'],
        loci=broadcast_values(
            params['phenome']['loci'],
            params['phenome']['n_phenotypes']
        ),
        h2_genetic=broadcast_values(
            params['phenome']['h2_genetic'],
            params['phenome']['n_phenotypes']
        ),
        h2_environmental=broadcast_values(
            params['phenome']['h2_environmental'],
            params['phenome']['n_phenotypes']
        ),
        h2_vertical=broadcast_values(
            params['phenome']['h2_vertical'],
            params['phenome']['n_phenotypes']
        ),
        genetic_cor=column_major(params['phenome']['genetic_cor']),
        environmental_cor=column_major(
            params['phenome']['environmental_cor']
        ),
    )

    b.mating(
        mating_type=params['mating']['mating_type'],
        n_iterations=params['mating']['n_iterations'],
        temp_init=params['mating']['temp_init'],
        temp_decay=params['mating']['temp_decay'],
        mate_cor=column_major(params['mating']['mate_cor']),
    )

    metric_specs = build_metrics(params['metrics']['metrics'])
    b.metrics(metric_specs=metric_specs)

    sim = b.build()
    sim.run()


class Simulation:
    def __init__(self) -> None:
        self._params: dict = {}

    def simulation(
        self,
        n_generations: int,
        n_individuals: int,
        output_dir: str | Path,
        random_seed: int | None = None,
    ) -> Self:
        if random_seed is None:
            random_seed = int(time_ns()) & 0xFFFFFFFFF

        path = Path(output_dir)
        path.mkdir(parents=True, exist_ok=True)
        output_dir = str(path.resolve())

        self._params['simulation'] = {
            'n_generations': n_generations,
            'n_individuals': n_individuals,
            'output_dir': output_dir,
            'random_seed': random_seed,
        }
        return self

    def genome(
        self,
        n_loci: int,
        locus_mafs: float | list[float],
        locus_recombination: float | list[float],
        locus_mutation: float | list[float],
    ) -> Self:
        self._params['genome'] = {
            'n_loci': n_loci,
            'locus_mafs': locus_mafs,
            'locus_recombination': locus_recombination,
            'locus_mutation': locus_mutation,
        }
        return self

    def phenome(
        self,
        n_phenotypes: int,
        names: list[str],
        loci: int | list[int],
        h2_genetic: float | list[float],
        h2_environmental: float | list[float],
        h2_vertical: float | list[float],
        genetic_cor: list[list[float]] | np.ndarray,
        environmental_cor: list[list[float]] | np.ndarray,
    ) -> Self:
        self._params['phenome'] = {
            'n_phenotypes': n_phenotypes,
            'names': names,
            'loci': loci,
            'h2_genetic': h2_genetic,
            'h2_environmental': h2_environmental,
            'h2_vertical': h2_vertical,
            'genetic_cor': genetic_cor,
            'environmental_cor': environmental_cor,
        }
        return self

    def mating(
        self,
        mating_type: MatingType,
        n_iterations: int,
        temp_init: float,
        temp_decay: float,
        mate_cor: list[list[float]] | np.ndarray,
    ) -> Self:
        self._params['mating'] = {
            'mating_type': mating_type,
            'n_iterations': n_iterations,
            'temp_init': temp_init,
            'temp_decay': temp_decay,
            'mate_cor': mate_cor,
        }
        return self

    def metrics(self, metrics: list[Metric]) -> Self:
        self._params['metrics'] = {'metrics': metrics}
        return self

    def run(self, n_replicates: int = 1, n_workers: int = 1) -> None:
        """Run one or multiple simulation replicates."""
        if not self._params:
            raise RuntimeError("Simulation parameters have not been set.")

        if n_replicates == 1:
            simulate_replication(self._params)
        else:
            with Pool(processes=n_workers) as pool:
                pool.starmap(
                    simulate_replication,
                    [
                        (self._params, rep_id)
                        for rep_id in range(n_replicates)
                    ],
                )

