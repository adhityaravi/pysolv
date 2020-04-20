"""Jacobi linear solver.

Jacobi solver is a type of stationary iteration scheme for solving a set of linear equations. The general linear
iteration scheme reads:

x(k+1) = x(k) + (1/C)(b - Ax(k))

where A is the coefficient matrix
      b is the RHS vector
      x is the solution vector
      C is the approximation of A that is easier to invert
      k is the iteration number

In the Jacobi scheme, C is approximated by only the diagonal values of A. That is,

A = L + D + U; C = D

Jacobi is a simple iterative scheme that converges only if the A matrix is strictly and irreducibly diagonally dominant.
Or, mathematically, Jacobi scheme converges only if the spectral radius of the matrix (1/D)(L+U) is < 1. pysolv
provides a function for checking the diagonal dominance of the A matrix.
"""

# import the necessary packages
from data import *


class JacobiSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
        """

        # data inheritance
        Data.__init__(self)

        # call the jacobi solver
        self.solve()

    def solve(self):
        """Jacobi solver.
        """

        print('wubba lubba dub dub')