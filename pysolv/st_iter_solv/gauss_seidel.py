"""Gauss-Seidel linear solver.

Gauss-Seidel is a type of stationary iteration scheme for solving a set of linear equations. The general linear
iteration scheme reads:

x(k+1) = x(k) + (1/C)(b - Ax(k))

where A is the coefficient matrix
      b is the RHS vector
      x is the solution vector
      C is the approximation of A that is easier to invert
      k is the iteration number

In the Gauss-Seidel scheme, C is approximated by combining the lower and upper triangular matrices of A. That is,

A = L + D + U; C = L + U

Initial guess for the solution vector, if not mentioned, is taken as an unity vector.

Properties:

1. The Gauss-Seidel scheme converges if either,

    i. The A matrix is strictly or irreducibly diagonally dominant. Or

    ii. The A matrix is symmetric positive-definite (symmetric --> A = transpose(A); positive definite --> all the
        eigen values of A are positive)

2. The convergence of the Gauss-Seidel scheme can be checked posteriori using the condition: spectral radius of
   (1/(L+U))(D) is < 1.

3. Conventional Gauss-Seidel scheme cannot handle non-square A matrices. An updated Gauss-Seidel scheme called
   generalized Gauss-Seidel is required.
"""

# import the necessary packages
from pysolv.data import *
from pysolv.f_lib import f_sorsolve
import numpy as np


class GaussSeidelSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
        """

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self._init_x0()

        # call the gauss seidel solver
        self._solve()

    def _init_x0(self):
        """Check if an initial value exists for the solution vector or initialize the solution vector"""

        try:
            self.x0
        except AttributeError:
            self.x0 = np.ones(self.n)

        # initialize the solution vectors
        self.x = np.empty(self.n)   # solution vector at the i'th iteration

    def _solve(self):
        """Serial python implementation of Gauss Seidel solver."""

        # currently cannot handle non square systems in Gauss-Seidel
        if self.m != self.n:
            msg = 'pysolv does not include a gauss-seidel solver for non-square systems yet'
            raise ValueError(msg)

        # Gauss-Seidel iterations
        # x_i(iter) = (1/a_i_i)(b_i - sum_(j=1)^(i-1)(a_i_j * x_j(iter)) - sum_(j=i+1)^(n)(A_i_j * x_j(iter-1)))
        else:
            # set relaxation parameter = 1
            omega = 1
            # dummy parameters
            h = 2
            adaptive_omega = 0
            omega_update_frequency = 1
            # call the sor fortran wrapper with relaxation parameter = 1
            f_sorsolve.init(self.A, self.b, self.x0, self.itermax, self.tol, omega, h, adaptive_omega,
                            omega_update_frequency, self.c1, self.c2, self.lambda1, self.lambda2, self.rho1)
            f_sorsolve.solve(self.x)

            # add the solution to the Data class
            Data.add_data('x', self.x)
