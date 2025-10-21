from typing import Union, List

import numpy as np

def broadcast_values(value: Union[float, List[float]], n_elems: int) -> List[float]:
    """
    Broadcast a single value or validate list length

    Parameters
    ----------
    value : float or list of float
        Single value to broadcast, or a list
    n_elems : int
        Expected length of output list

    Returns
    -------
    list of float
        List of length `n_elems`
    """
    if isinstance(value, (int, float)):
        return [float(value)] * n_elems
    elif isinstance(value, list):
        if len(value) != n_elems:
            raise ValueError(f'expected {n_elems} values; got {len(value)}')
        return value 
    else:
        raise TypeError(f'expected float or list, got {type(value)}')

def column_major(matrix: Union[List[List[float]], np.ndarray]) -> List[float]:
    """
    Convert a row-major order matrix to column-major

    Parameters
    ----------
    matrix : array_like
        2D major in row-major order

    Returns
    -------
    list of float
        Flattened matrix in column-major (FORTRAN) order
    """
    arr = np.asarray(matrix, dtype=float)
    return arr.flatten(order='F').tolist()

def check_cor(matrix: Union[List[List[float]], np.ndarray]) -> bool:
    """
    Check whether a matrix is a valid correlation matrix

    Parameters
    ----------
    matrix: array_like
        2D matrix in row-major order to check
    """
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.all(eigenvalues >= 0).item()
