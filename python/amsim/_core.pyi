"""Type stubs for C++ _core module"""


class MatingType:
    RANDOM: MatingType
    ASSORTATIVE: MatingType


class ComponentType:
    GENETIC: ComponentType
    ENVIRONMENTAL: ComponentType
    VERTICAL: ComponentType
    TOTAL: ComponentType


class _MetricSpec:
    def __init__(self) -> None: ...


# register these in backend...
# def _loc_maf(n_loci: int) -> MetricSpec: ...
# def _loc_mean(n_loci: int) -> MetricSpec: ...
# def _loc_var(n_loci: int) -> MetricSpec: ...
def _pheno_h2() -> _MetricSpec: ...
def _pheno_comp_mean(component_type: ComponentType) -> _MetricSpec: ...
def _pheno_comp_var(component_type: ComponentType) -> _MetricSpec: ...
def _pheno_comp_cor(component_type: ComponentType) -> _MetricSpec: ...
def _pheno_comp_xcor(component_type: ComponentType) -> _MetricSpec: ...
def _pheno_latent_h2() -> _MetricSpec: ...
def _pheno_latent_comp_mean(component_type: ComponentType) -> _MetricSpec: ...
def _pheno_latent_comp_var(component_type: ComponentType) -> _MetricSpec: ...
def _pheno_latent_comp_cor(component_type: ComponentType) -> _MetricSpec: ...
def _pheno_latent_comp_xcor(component_type: ComponentType) -> _MetricSpec: ...

class _SimulationBuilder:
    def __init__(self) -> None: ...

    def simulation(
        self,
        n_generations: int,
        n_individuals: int,
        output_dir: str,
        random_seed: int
    ) -> _SimulationBuilder: ...

    def genome(
        self,
        n_loci: int,
        locus_mafs: list[float],
        locus_recombination: list[float],
        locus_mutation: list[float]
    ) -> _SimulationBuilder: ...

    def phenome(
        self,
        n_phenotypes: int,
        names: list[str],
        loci: list[int],
        h2_genetic: list[float],
        h2_environmental: list[float],
        h2_vertical: list[float],
        genetic_cor: list[float],
        environmental_cor: list[float]
    ) -> _SimulationBuilder: ...

    def mating(
        self,
        mating_type: MatingType,
        n_iterations: int | None,
        temp_init: float | None,
        temp_decay: float | None,
        mate_cor: list[float]
    ) -> _SimulationBuilder: ...

    def metrics(
        self,
        metric_specs: list[_MetricSpec]
    ) -> _SimulationBuilder: ...

    def build(self) -> _Simulation: ...

class _Simulation:
    def run(self) -> None: ...

