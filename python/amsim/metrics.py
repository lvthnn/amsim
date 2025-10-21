from amsim._core import ComponentType

from amsim._core import (
    _Metric,
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

def loc_maf(n_loci: int) -> _Metric:
    """
    Create locus MAF metric.

    Parameters
    ----------
    n_loci : int
        Number of simulated loci
    """
    return _loc_maf(n_loci)

def loc_mean(n_loci: int) -> _Metric:
    """
    Create locus mean metric.

    Parameters
    ----------
    n_loci : int
        Number of simulated loci
    """
    return _loc_mean(n_loci)

def loc_var(n_loci: int) -> _Metric:
    """
    Create locus variance metric.

    Parameters
    ----------
    n_loci : int
        Number of simulated loci
    """
    return _loc_var(n_loci)

def pheno_h2(n_phenotypes: int) -> _Metric:
    """
    Create phenotype heritability metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    """
    return _pheno_h2(n_phenotypes)

def pheno_comp_mean(n_phenotypes: int, component_type: ComponentType) -> _Metric:
    """
    Create phenotype component mean metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return _pheno_comp_mean(n_phenotypes, component_type)

def pheno_comp_var(n_phenotypes: int, component_type: ComponentType) -> _Metric:
    """
    Create phenotype component variance metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return _pheno_comp_var(n_phenotypes, component_type)

def pheno_comp_cor(n_phenotypes: int, component_type: ComponentType) -> _Metric:
    """
    Create phenotype component correlation metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return _pheno_comp_cor(n_phenotypes, component_type)

def pheno_comp_xcor(n_phenotypes: int, component_type: ComponentType) -> _Metric:
    """
    Create phenotype component cross-correlation metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return _pheno_comp_xcor(n_phenotypes, component_type)

def pheno_latent_h2(n_phenotypes: int) -> _Metric:
    """
    Create latent phenotype heritability metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    """
    return _pheno_latent_h2(n_phenotypes)

def pheno_latent_comp_cor(n_phenotypes: int, component_type: ComponentType) -> _Metric:
    """
    Create latent phenotype component correlation metric.

    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return _pheno_latent_comp_cor(n_phenotypes, component_type)

def pheno_latent_comp_xcor(n_phenotypes: int, component_type: ComponentType) -> _Metric:
    """
    Create latent phenotype component cross-correlation metric.
    
    Parameters
    ----------
    n_phenotypes : int
        Number of simulated phenotypes
    component_type : ComponentType
        Phenotype component type
    """
    return _pheno_latent_comp_xcor(n_phenotypes, component_type)
