from typing import TypeVar

import numpy as np

T = TypeVar('T', int, float)

def broadcast_values(values: T | list[T], n_elems: int) -> list[T]:
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
    if isinstance(values, list):
        if len(values) != n_elems:
            raise ValueError(f'expected {n_elems} values; got {len(values)}')
        return values 
    else:
        return [values] * n_elems

def column_major(matrix: list[list[T]] | np.ndarray) -> list[T]:
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
    arr = np.asarray(matrix)
    return arr.flatten(order='F').tolist()

def check_cor(matrix: list[list[T]] | np.ndarray) -> bool:
    """
    Check whether a matrix is a valid correlation matrix

    Parameters
    ----------
    matrix: array_like
        2D matrix in row-major order to check
    """
    eigenvalues = np.linalg.eigvalsh(matrix)
    matrix_np = np.asarray(matrix)
    val = np.all(matrix_np >= 0).item() and np.all(matrix_np <= 1).item()
    psd = np.all(eigenvalues >= 0).item()
    return val and psd
