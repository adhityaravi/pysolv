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
from data import *
import numpy as np
import time as ti


class GaussSeidelSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
           ToDO: Deprecate Gauss-Seidel _solve() and rather call SOR with omega=1
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
        self.x_old = self.x0  # solution vector at the i-1'th iteration

    def _solve(self):
        """Serial python implementation of Gauss Seidel solver."""

        # currently cannot handle non square systems in Gauss-seidel
        # ToDo: Implement generalized Gauss-seidel for non-square systems based on M.Saha's paper
        if self.m != self.n:
            msg = 'pysolv does not include a gauss-seidel solver for non-square systems yet'
            raise ValueError(msg)

        # Gauss-Seidel iterations
        # x_i(iter) = (1/a_i_i)(b_i - sum_(j=1)^(i-1)(a_i_j * x_j(iter)) - sum_(j=i+1)^(n)(A_i_j * x_j(iter-1)))
        else:
            # time the solver
            start_time = ti.time()

            # Initialize the iteration counter and solution vector
            iter = 0

            # continue the iteration till the iteration counter reaches the maximum count if convergence is not obtained
            while iter < self.ITERMAX:
                for i in range(0, self.n):
                    s1 = 0
                    s2 = 0

                    for j in range(0, i):
                        s1 = s1 + (self.A[i, j]*self.x[j])

                    for j in range(i+1, self.n):
                        s2 = s2 + (self.A[i, j]*self.x_old[j])

                    self.x[i] = (1/self.A[i, i])*(self.b[i] - s1 - s2)

                # stopping criteria: (||x(iter) - x(iter-1)|| / ||x(iter)||) < TOL
                # ToDo: Implement a pysolv native function to compute the residual
                res = np.linalg.norm((self.x - self.x_old)) / np.linalg.norm(self.x_old)
                if res < self.TOL:
                    break
                # Update the iteration counter and x(iter-1)
                self.x_old = self.x.copy()
                iter += 1

            # time taken for convergence
            time_taken = ti.time() - start_time

            # add the solution to the Data class
            Data.add_data('x', self.x)
            Data.add_data('time_taken', time_taken)
            Data.add_data('iterations', iter - 1)
            Data.add_data('residual', res)