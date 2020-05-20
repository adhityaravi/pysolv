"""Pre-Conditioners for Krylov subspace based linear solvers.

ToDo: Current approaches are naive and too expensive. Improve them using approaches from Yousef Saad.
"""

# import the necessary packages
import numpy as np


def jacobi(A):
    """Jacobi pre-conditioner.

    Jacobi pre-conditioner is one of the simplest pre-conditioner used for improving the convergence of different
    non-stationary iterative solvers by improving the condition number of the coefficient matrix.

    Jacobi pre-conditioner is a diagonal matrix whose diagonal elements are same as that of the input coefficient
    matrix.

    M = diag(A)

    Parameters:
    ----------
        A (numpy array): coefficient matrix

    Returns:
    -------
        M (numpy array): Jacobi pre-conditioner
    """

    return np.diag(np.diag(A))


def gauss_seidel(A):
    """Gauss-Seidel pre-conditioner.

    Gauss-Seidel or the SOR (with omega = 1) pre-conditioner is a pre-conditioner used for improving the convergence of
    different non-stationary iterative solvers by improving the condition number of the coefficient matrix.

    The omega value for the SOR pre-conditioner is chosen as 1 (--> Gauss-Seidel precon) based on Yousef Saad's book.

    The Gauss-Seidel pre-conditioner is built using the lower triangular matrix of the coefficient matrix

    A = L + D + U
    M = L + D

    Parameters:
    ----------
        A (numpy array): coefficient matrix

    Returns:
    -------
        M (numpy array): Gauss-Seidel pre-conditioner
    """

    return np.tril(A)
