"""Implements pickleable metric specification for parallel simulations"""
from amsim._core import(
    _MetricSpec,
    _pheno_h2,
    _pheno_comp_mean,
    _pheno_comp_var,
    _pheno_comp_cor,
    _pheno_comp_xcor,
    _pheno_latent_h2,
    _pheno_latent_comp_mean,
    _pheno_latent_comp_var,
    _pheno_latent_comp_cor,
    _pheno_latent_comp_xcor
)

from amsim._core import ComponentType 

_METRICS = {
    "pheno_h2": _pheno_h2,
    "pheno_comp_mean": _pheno_comp_mean,
    "pheno_comp_var": _pheno_comp_var,
    "pheno_comp_cor": _pheno_comp_cor,
    "pheno_comp_xcor": _pheno_comp_xcor,
    "pheno_latent_h2": _pheno_latent_h2,
    "pheno_latent_comp_mean": _pheno_latent_comp_mean,
    "pheno_latent_comp_var": _pheno_latent_comp_var,
    "pheno_latent_comp_cor": _pheno_latent_comp_cor,
    "pheno_latent_comp_xcor": _pheno_latent_comp_xcor
}

class Metric:
    _name: str
    _args: tuple | None

    def __init__(self, name: str, args: tuple | None = None) -> None:
        self._name, self._args = name, args

    def build(self) -> _MetricSpec:
        if self._args is None:
            return _METRICS[self._name]()
        else:
            return _METRICS[self._name](*self._args)

def build_metrics(metric_specs: list[Metric]) -> list[_MetricSpec]:
    return [metric_spec.build() for metric_spec in metric_specs]

def pheno_h2():
    return Metric("pheno_h2")

def pheno_comp_mean(component_type: ComponentType):
    return Metric("pheno_comp_mean", (component_type, ))

def pheno_comp_var(component_type: ComponentType):
    return Metric("pheno_comp_var", (component_type, ))

def pheno_comp_cor(component_type: ComponentType):
    return Metric("pheno_comp_cor", (component_type, ))

def pheno_comp_xcor(component_type: ComponentType):
    return Metric("pheno_comp_xcor", (component_type, ))

def pheno_latent_h2():
    return Metric("pheno_latent_h2")

def pheno_latent_comp_mean(component_type: ComponentType):
    return Metric("pheno_latent_comp_mean", (component_type, ))

def pheno_latent_comp_var(component_type: ComponentType):
    return Metric("pheno_latent_comp_var", (component_type, ))

def pheno_latent_comp_cor(component_type: ComponentType):
    return Metric("pheno_latent_comp_cor", (component_type, ))

def pheno_latent_comp_xcor(component_type: ComponentType):
    return Metric("pheno_latent_comp_xcor", (component_type, ))

