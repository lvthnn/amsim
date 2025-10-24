"""Implements user interface to amsimcpp backend"""
from time import time_ns
from pathlib import Path
from typing import Self
from multiprocessing import Pool
from configparser import ConfigParser

import numpy as np

from amsim._core import _SimulationBuilder, _Simulation
from amsim._core import MatingType
from amsim.utils import broadcast_values, column_major
from amsim.metrics import Metric, build_metrics

def simulate_replication(params: dict, rep_id : int | None = None) -> None:
    if rep_id is not None:
        output_dir = Path(params['simulation']['output_dir'])
        rep_dir = output_dir / f'rep_{rep_id:03d}'
        rep_dir.mkdir(parents=True, exist_ok=True)
        params['simulation']['output_dir'] = str(rep_dir.resolve())
        params['simulation']['random_seed'] = params['simulation']['random_seed'] + rep_id

    builder = SimulationBuilder().simulation(
        n_generations=params['simulation']['n_generations'],
        n_individuals=params['simulation']['n_individuals'],
        output_dir=params['simulation']['output_dir'],
        random_seed=params['simulation']['random_seed']
    ).genome(
        n_loci=params['genome']['n_loci'],
        locus_mafs=params['genome']['locus_mafs'],
        locus_recombination=params['genome']['locus_recombination'],
        locus_mutation=params['genome']['locus_mutation']
    ).phenome(
        n_phenotypes=params['phenome']['n_phenotypes'],
        names=params['phenome']['names'],
        loci=params['phenome']['loci'],
        h2_genetic=params['phenome']['h2_genetic'],
        h2_environmental=params['phenome']['h2_environmental'],
        h2_vertical=params['phenome']['h2_vertical'],
        genetic_cor=params['phenome']['genetic_cor'],
        environmental_cor=params['phenome']['environmental_cor']
    ).mating(
        mating_type=params['mating']['mating_type'],
        n_iterations=params['mating']['n_iterations'],
        temp_init=params['mating']['temp_init'],
        temp_decay=params['mating']['temp_decay'],
        mate_cor=params['mating']['mate_cor']
    ).metrics(
        metrics=params['metrics']['metrics']
    )

    simulation = builder.build()
    simulation.run()

    pass

class SimulationBuilder:
    def __init__(self) -> None:
        self._builder = _SimulationBuilder()

    def simulation(
        self,
        n_generations : int,
        n_individuals : int,
        output_dir : str | Path,
        random_seed : int | None
    ) -> Self:
        if random_seed is None:
            random_seed = int(time_ns()) & 0xFFFFFFFFF

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        if isinstance(output_dir, Path):
            output_dir = str(output_dir.resolve())

        self._builder.simulation(
            n_generations=n_generations,
            n_individuals=n_individuals,
            output_dir=output_dir,
            random_seed=random_seed
        )
        return self

    def genome(
        self,
        n_loci: int,
        locus_mafs: float | list[float],
        locus_recombination: float | list[float],
        locus_mutation: float | list[float]
    ) -> Self:
        self._builder.genome(
            n_loci=n_loci,
            locus_mafs=broadcast_values(locus_mafs, n_loci),
            locus_recombination=broadcast_values(locus_recombination, n_loci),
            locus_mutation=broadcast_values(locus_mutation, n_loci)
        )
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
        environmental_cor: list[list[float]] | np.ndarray
    ) -> Self:
        self._builder.phenome(
            n_phenotypes=n_phenotypes,
            names=names,
            loci=broadcast_values(loci, n_phenotypes),
            h2_genetic=broadcast_values(h2_genetic, n_phenotypes),
            h2_environmental=broadcast_values(h2_environmental, n_phenotypes),
            h2_vertical=broadcast_values(h2_vertical, n_phenotypes),
            genetic_cor=column_major(genetic_cor),
            environmental_cor=column_major(environmental_cor)
        )
        return self

    def mating(
        self,
        mating_type: MatingType,
        n_iterations: int,
        temp_init: float,
        temp_decay: float,
        mate_cor: list[list[float]] | np.ndarray
    ) -> Self:
        self._builder.mating(
            mating_type=mating_type,
            n_iterations=n_iterations,
            temp_init=temp_init,
            temp_decay=temp_decay,
            mate_cor=column_major(mate_cor)
        )
        return self

    def metrics(self, metrics: list[Metric]) -> Self:
        _metric_specs = build_metrics(metrics)
        self._builder.metrics(metric_specs=_metric_specs)
        return self

    def build(self) -> _Simulation:
        return self._builder.build()


class Simulation:
    def __init__(self) -> None:
        self._params = dict()

    def simulation(
        self,
        n_generations : int,
        n_individuals : int,
        output_dir : str | Path,
        random_seed : int | None
    ) -> Self:
        if random_seed is None:
            random_seed = int(time_ns()) & 0xFFFFFFFFF

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        if isinstance(output_dir, Path):
            output_dir = str(output_dir.resolve())

        self._params['simulation'] = {
            'n_generations': n_generations,
            'n_individuals': n_individuals,
            'output_dir': output_dir,
            'random_seed': random_seed
        }
        return self

    def genome(
        self,
        n_loci: int,
        locus_mafs: float | list[float],
        locus_recombination: float | list[float],
        locus_mutation: float | list[float]
    ) -> Self:
        self._params['genome'] = {
            'n_loci': n_loci,
            'locus_mafs': broadcast_values(locus_mafs, n_loci),
            'locus_recombination': broadcast_values(locus_recombination, n_loci),
            'locus_mutation': broadcast_values(locus_mutation, n_loci)
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
        environmental_cor: list[list[float]] | np.ndarray
    ) -> Self:
        self._params['phenome'] = {
            'n_phenotypes': n_phenotypes,
            'names': names,
            'loci': broadcast_values(loci, n_phenotypes),
            'h2_genetic': broadcast_values(h2_genetic, n_phenotypes),
            'h2_environmental': broadcast_values(h2_environmental, n_phenotypes),
            'h2_vertical': broadcast_values(h2_vertical, n_phenotypes),
            'genetic_cor': column_major(genetic_cor),
            'environmental_cor': column_major(environmental_cor)
        }    
        return self

    def mating(
        self,
        mating_type: MatingType,
        n_iterations: int,
        temp_init: float,
        temp_decay: float,
        mate_cor: list[list[float]] | np.ndarray
    ) -> Self:
        self._params['mating'] = {
            'mating_type': mating_type,
            'n_iterations': n_iterations,
            'temp_init': temp_init,
            'temp_decay': temp_decay,
            'mate_cor': column_major(mate_cor)
        }
        return self

    def metrics(self, metrics: list[Metric]) -> Self:
        self._params['metrics'] = {
            'metrics': metrics 
        }
        return self

    def run(self, n_replicates: int = 1, n_threads: int = 1):
        if n_replicates == 1:
            simulate_replication(self._params)
        else:
            with Pool(processes=n_threads) as pool:
                pool.starmap(
                    simulate_replication,
                    [(self._params, rep_id) for rep_id in range(n_replicates)]
                )

