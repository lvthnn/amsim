"""Implements user interface to amsimcpp backend"""
from time import time_ns
from pathlib import Path
from typing import Self
from multiprocessing import Pool
from configparser import ConfigParser

import numpy as np

from amsim._core import _SimulationBuilder, _Simulation, MatingType
from amsim.metrics import MetricSpec, build_metrics

from amsim.utils import broadcast_values, column_major


def _run_single_simulation(args: tuple[dict, int]):
    params, rep_id = args

    from amsim import SimulationBuilder

    output_rep = Path(params['simulation']['output_dir']) / f'rep_{(rep_id + 1):03d}'
    output_rep.mkdir(parents=True, exist_ok=True)

    builder = SimulationBuilder()

    builder.simulation(
        n_generations=params['simulation']['n_generations'],
        n_individuals=params['simulation']['n_individuals'],
        output_dir=output_rep,
        random_seed=None
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
        metric_specs=params['metrics']['metric_specs']
    )

    simulation = builder.build()
    simulation.run()

    return output_rep


class Simulation:
    def __init__(self, params: dict, simulation: _Simulation) -> None:
        self._params = params
        self._simulation = simulation
        self._results = None

    def run(self, n_replicates: int = 1, n_threads: int = 1) -> None:
        params = self._params.copy()
        params['metrics'] = {
            spec._name: spec._args for spec in params['metrics']['metric_specs']
        }
        config = ConfigParser()
        config_path = Path(f'{params['simulation']['output_dir']}/params.cfg')
        config.read_dict(params)

        with open(config_path, 'w') as config_file:
            config.write(config_file)
        
        if n_replicates == 1:
            self._simulation.run()
        else:
            self._run_replicates(n_replicates, n_threads)

    def _run_replicates(self, n_replicates: int, n_threads: int) -> None:
        output_dirs = []
        args = [(self._params, rep) for rep in range(n_replicates)]

        with Pool(processes=n_threads) as pool:
            output_dirs = pool.map(_run_single_simulation, args)

        self._results = output_dirs


class SimulationBuilder:
    def __init__(self) -> None:
        self._params = dict()
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

        self._params['simulation'] = {
            'n_generations': n_generations,
            'n_individuals': n_individuals,
            'output_dir': output_dir,
            'random_seed': random_seed
        }

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
        self._params['genome'] = {
            'n_loci': n_loci,
            'locus_mafs': broadcast_values(locus_mafs, n_loci),
            'locus_recombination': broadcast_values(locus_recombination, n_loci),
            'locus_mutation': broadcast_values(locus_mutation, n_loci)
        }

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
        self._params['mating'] = {
            'mating_type': mating_type,
            'n_iterations': n_iterations,
            'temp_init': temp_init,
            'temp_decay': temp_decay,
            'mate_cor': column_major(mate_cor)
        }

        self._builder.mating(
            mating_type=mating_type,
            n_iterations=n_iterations,
            temp_init=temp_init,
            temp_decay=temp_decay,
            mate_cor=column_major(mate_cor)
        )
        return self

    def metrics(self, metric_specs: list[MetricSpec]) -> Self:
        self._params['metrics'] = { 'metric_specs': metric_specs }
        metrics = build_metrics(metric_specs)
        self._builder.metrics(metrics=metrics)
        return self

    def build(self) -> Simulation:
        return Simulation(self._params, self._builder.build())
