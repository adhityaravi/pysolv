"""Pysolv solve.

solve() is used to solve a set of linear equations.

solve() uses various stationary iterative and Krylov subspace based schemes for solving the linear set of equations.
"""

# import the necessary packages
from pysolv.data import Data
from pysolv.tools import *


def solve(A, b, solver='jacobi', **kwargs):
    """Solve a linear system of equations

    Parameters:
         A (numpy array): coefficient matrix
         b (numpy array): RHS vector
         solver (string): Type of the scheme to be used to solve the equation, by default uses Jacobi

    Returns:
        x (numpy array): solution to the linear system
        res (float): residual with which the linear solver converged
        it (int): number of iterations at which the linear solver converged
        time_taken (float): time taken by the linear solver to converge in [s]
    """

    # add the coefficient matrix, RHS vector to the Data class
    Data.add_data('A', A)
    Data.add_data('b', b)

    # check if the solver scheme is available in pysolv and add it to the Data class
    check_solver(solver)

    # add the rest of the keyword arguments, if any, to the Data class
    add_kwargs(kwargs)

    # extract different properties of the linear system and store it in the Data class
    ext_lsprops()

    # call the appropriate solver to solve the linear system
    call_solver()

    # fetch the solution from the Data class
    x = Data.x

    # flush the Data class clean
    Data._Data__reset_data()

    # return the solution
    return x
