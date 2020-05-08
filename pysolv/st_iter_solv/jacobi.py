"""Jacobi linear solver.

Jacobi is a type of stationary iteration scheme for solving a set of linear equations. The general linear
iteration scheme reads:

x(k+1) = x(k) + (1/C)(b - Ax(k))

where A is the coefficient matrix
      b is the RHS vector
      x is the solution vector
      C is the approximation of A that is easier to invert
      k is the iteration number

In the Jacobi scheme, C is approximated by only the diagonal values of A. That is,

A = L + D + U; C = D

Initial guess for the solution vector, if not mentioned, is taken as an unity vector.

Properties:

1. Jacobi is a very simple iterative scheme that works if and only if all the diagonal elements of the A matrix
   are non zero.

2. Converges only if the A matrix is strictly and irreducibly diagonally dominant. Or, mathematically,
   Jacobi scheme converges only if the spectral radius of the matrix (1/D)(L+U) is < 1.

3. Can become memory expensive for large systems as duplication of the solution vector x(k) is necessary.

4. Conventional Jacobi scheme cannot handle non-square A matrices. An updated Jacobi scheme called generalized
   Jacobi is required.
"""

# import the necessary packages
from pysolv.data import *
from pysolv.f_lib import f_jacobisolve
import numpy as np
import time as ti


class JacobiSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
        """

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self._init_x0()

        # check for relaxation parameter
        self._check_omega()

        # call the jacobi solver
        self._solve()

    def _init_x0(self):
        """Check if an initial value exists for the solution vector or initialize the solution vector
        """

        try:
            self.x0
        except AttributeError:
            self.x0 = np.ones(self.n)

        # initialize the solution vectors
        self.x = np.empty(self.n)  # solution at the i'th iteration

    def _check_omega(self):
        """Check if the user has prescribed a value for the relaxation parameter. If not, prescribe a value of 1 (Jacobi
        iterations are un-damped by default)
        """

        try:
            self.omega
        except AttributeError:
            self.omega = 1

    def _solve(self):
        """Serial python implementation of conventional Jacobi solver.
        """

        # currently cannot handle non square systems in Jacobi
        # ToDo: Implement generalized Jacobi for non-square systems based on M.Saha's paper
        if self.m != self.n:
            msg = 'pysolv does not include a jacobi solver for non-square systems yet'
            raise ValueError(msg)

        # Jacobi iteration:
        # x_i(iter) = (1/a_i_i)(b_i - sum(a_i_j * x_j(iter-1))) where, i = 1 to n, j = 1 to n & j != i
        else:
            # call the fortran jacobi solver
            f_jacobisolve.solve(self.A, self.b, self.x0, self.x, self.itermax, self.tol, self.omega)

            # add the solution to the Data class
            Data.add_data('x', self.x)
