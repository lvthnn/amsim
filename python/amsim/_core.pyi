"""Type stubs for C++ _core module"""

class SimulationBuilder:
    def __init__(self) -> None: ...

    def simulation(
        self,
        n_generations: int,
        n_individuals: int,
        output_dir:    str,
        random_seed:   int
    ) -> SimulationBuilder: ...

    def genome(
        self,
        n_loci:              int,
        locus_mafs:          list[float],
        locus_recombination: list[float],
        locus_mutate:        list[float]
    ) -> SimulationBuilder: ...

    def phenome(
        self,
        n_phenotypes:      int,
        names:             list[str],
        loci:              list[int],
        h2_genetic:        list[float],
        h2_environmental:  list[float],
        h2_vertical:       list[float],
        genetic_cor:       list[float],
        environmental_cor: list[float]
    ) -> SimulationBuilder: ...

    def mating(
        self,
        mating_type:  MatingType,
        n_iterations: int | None,
        temp_init:    float | None,
        temp_decay:   float | None,
        mate_cor:     list[float]
    ) -> SimulationBuilder: ...

    def metric(
        self,
        metrics: list[Metric]
    ) -> SimulationBuilder: ...

    def build(self) -> Simulation: ...

class Simulation:
    def run(self) -> None: ...

class MatingType:
    RANDOM: MatingType
    ASSORTATIVE: MatingType

class Metric:
    def __init__(self) -> None: ...
