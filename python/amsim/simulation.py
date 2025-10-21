from pathlib import Path
from time import time_ns
from typing import Self, List

from amsim._core import _SimulationBuilder, _Simulation, MatingType
from amsim._core import _Metric

from amsim.utils import broadcast_values

class SimulationBuilder:
    def __init__(self):
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

        self._builder.simulation(n_generations=n_generations,
                                 n_individuals=n_individuals,
                                 output_dir=str(output_dir),
                                 random_seed=random_seed)
        return self

    def genome(self, n_loci: int, locus_mafs: float | List[float],
               locus_recombination: float | List[float],
               locus_mutation: float | List[float]) -> Self:
        self._builder.genome(n_loci=n_loci,
                             locus_mafs=broadcast_values(locus_mafs, n_loci),
                             locus_recombination=broadcast_values(locus_recombination, n_loci),
                             locus_mutation=broadcast_values(locus_mutation, n_loci))
        return self

    def phenome(self, n_phenotypes: int, names: list[str], loci: list[int],
                h2_genetic: list[float], h2_environmental: list[float],
                h2_vertical: list[float], genetic_cor: list[float],
                environmental_cor: list[float]) -> Self:
        self._builder.phenome(n_phenotypes=n_phenotypes,
                              names=names,
                              loci=loci,
                              h2_genetic=h2_genetic,
                              h2_environmental=h2_environmental,
                              h2_vertical=h2_vertical,
                              genetic_cor=genetic_cor,
                              environmental_cor=environmental_cor)
        return self

    def mating(self, mating_type: MatingType, n_iterations: int,
               temp_init: float, temp_decay: float, mate_cor: list[float]) -> Self:
        self._builder.mating(mating_type=mating_type,
                             n_iterations=n_iterations,
                             temp_init=temp_init,
                             temp_decay=temp_decay,
                             mate_cor=mate_cor)
        return self


    def metric(self, metrics: list[_Metric]) -> Self:
        self._builder.metric(metrics=metrics)
        return self

    def build(self) -> _Simulation:
        return self._builder.build()

class Simulation:
    def __init__(self, config):
        raise NotImplementedError('config initialisation not supported yet')
