from amsim._core import ComponentType

from amsim._core import (
    _loc_maf,
    _loc_mean,
    _loc_var,
    _pheno_h2,
    _pheno_comp_mean,
    _pheno_comp_var,
    _pheno_comp_cor,
    _pheno_comp_xcor,
    _pheno_latent_h2,
    _pheno_latent_comp_cor,
    _pheno_latent_comp_xcor
) 

_METRIC_REGISTRY = {
    'loc_maf': _loc_maf,
    'loc_mean': _loc_mean,
    'loc_var': _loc_var,
    'pheno_h2': _pheno_h2,
    'pheno_comp_mean': _pheno_comp_mean,
    'pheno_comp_var': _pheno_comp_var,
    'pheno_comp_cor': _pheno_comp_cor,
    'pheno_comp_xcor': _pheno_comp_xcor,
    'pheno_latent_h2': _pheno_latent_h2,
    'pheno_latent_comp_cor': _pheno_latent_comp_cor,
    'pheno_latent_comp_xcor': _pheno_latent_comp_xcor
}

class MetricSpec:
    def __init__(self, name, args):
        self._name, self._args = name, args

def build_metrics(specs: list[MetricSpec]):
    return [_METRIC_REGISTRY[spec._name](*spec._args) for spec in specs]

def loc_maf(n_loci: int) -> MetricSpec:
    """
    Create locus MAF metric.

    Parameters
    ----------
    n_loci : int
        Number of simulated loci
    """
    return MetricSpec('loc_maf', (n_loci,))

def loc_mean(n_loci: int) -> MetricSpec:
    """
    Create locus mean metric.

    Parameters
    ----------
    n_loci : int
        Number of simulated loci
    """
    return MetricSpec('loc_mean', (n_loci,))

def loc_var(n_loci: int) -> MetricSpec:
    """
    Create locus variance metric.

    Parameters
    ----------
    n_loci : int
        Number of simulated loci
    """
    return MetricSpec('loc_var', (n_loci,))

def pheno_h2(n_phenotypes: int) -> MetricSpec:
    """
    Create phenotype heritability metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    """
    return MetricSpec("pheno_h2", (n_phenotypes,))

def pheno_comp_mean(n_phenotypes: int, component_type: ComponentType) -> MetricSpec:
    """
    Create phenotype component mean metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return MetricSpec("pheno_comp_mean", (n_phenotypes, component_type))

def pheno_comp_var(n_phenotypes: int, component_type: ComponentType) -> MetricSpec:
    """
    Create phenotype component variance metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return MetricSpec("pheno_comp_var", (n_phenotypes, component_type))

def pheno_comp_cor(n_phenotypes: int, component_type: ComponentType) -> MetricSpec:
    """
    Create phenotype component correlation metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return MetricSpec("pheno_comp_cor", (n_phenotypes, component_type))

def pheno_comp_xcor(n_phenotypes: int, component_type: ComponentType) -> MetricSpec:
    """
    Create phenotype component cross-correlation metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return MetricSpec("pheno_comp_xcor", (n_phenotypes, component_type))

def pheno_latent_h2(n_phenotypes: int) -> MetricSpec:
    """
    Create latent phenotype heritability metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    """
    return MetricSpec("pheno_latent_h2", (n_phenotypes,))

def pheno_latent_comp_cor(n_phenotypes: int, component_type: ComponentType) -> MetricSpec:
    """
    Create latent phenotype component correlation metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return MetricSpec("pheno_latent_comp_cor", (n_phenotypes, component_type))

def pheno_latent_comp_xcor(n_phenotypes: int, component_type: ComponentType) -> MetricSpec:
    """
    Create latent phenotype component cross-correlation metric.
    
    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return MetricSpec("pheno_latent_comp_xcor", (n_phenotypes, component_type))
