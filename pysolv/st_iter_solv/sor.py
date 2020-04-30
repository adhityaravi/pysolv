"""Successive over-relaxation linear solver.

Successive over-relaxation is a type of stationary iteration scheme for solving a set of linear equations.
The general linear iteration scheme reads:

x(k+1) = x(k) + (1/C)(b - Ax(k))

where A is the coefficient matrix
      b is the RHS vector
      x is the solution vector
      C is the approximation of A that is easier to invert
      k is the iteration number

In the Successive over-relaxation scheme, C is approximated by combining the lower and upper triangular matrices of A.
Then, a relaxation parameter omega is introduced for C. That is,

A = L + D + U; C = (1/omega) * (L + U)

Initial guess for the solution vector, if not mentioned, is taken as an unity vector.

Properties:
1. SOR is aimed at improving the convergence of the Gauss-Seidel solver by introducing the relaxation factor.

2. An apriori convergence criteria for the SOR solver is

    i. The A matrix is strictly or irreducibly diagonally dominant. Or

    ii. The A matrix is symmetric positive-definite (symmetric --> A = transpose(A); positive definite --> all the
        eigen values of A are positive)

3. A posteriori convergence criteria for the SOR scheme is: spectral radius of (1 - (1/omega*D + L)*A) is < 1.
"""

# import the necessary packages
from data import *
import numpy as np
import time as ti


class SORSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
        """

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self._init_x0()

        # check for relaxation parameter
        self._check_omega()

        # call the gauss seidel solver
        self._solve()

    def _init_x0(self):
        """Check if an initial value exists for the solution vector or initialize the solution vector
        """

        try:
            self.x0
        except AttributeError:
            self.x0 = np.ones(self.n)

        # initialize the solution vectors
        self.x = np.empty(self.n)   # solution vector at the i'th iteration
        self.x_old = self.x0  # solution vector at the i-1'th iteration

    def _check_omega(self):
        """Check if the user has prescribed a value for the relaxation parameter. If not, prescribe a value of 1.7 based
           on Dr. Rosics's lecture (Institute of Scientific Computing, TU Braunschweig):
           [https://www.tu-braunschweig.de/en/wire/teaching/previous-terms/winter-2016-17]
        """

        try:
            self.omega
        except AttributeError:
            self.omega = 1.7

    def _solve(self):
        """Serial python implementation of successive over-relaxation solver based on Dr. Rosic's lectures (Institute of
           Scientific Computing, TU Braunschweig):
           [https://www.tu-braunschweig.de/en/wire/teaching/previous-terms/winter-2016-17]
        """

        # currently the sor solver of pysolv cannot handle non-square systems
        if self.m != self.n:
            msg = 'pysolv does not include a sor solver for non-square systems yet'
            raise ValueError(msg)

        # SOR iteration
        # x_i(iter) = (1 - omega) * x_i(iter-1) +
        #             (omega/a_i_i)(b_i - sum_(j=1)^(i-1)(a_i_j * x_j(iter)) - sum_(j=i+1)^(n)(A_i_j * x_j(iter-1)))
        else:
            # time the solver
            start_time = ti.time()

            # initialize iteration counter
            iter = 0

            # continue the iteration till the iteration counter reaches the maximum count if convergence is not
            # obtained
            while iter < self.ITERMAX:
                for i in range(0, self.n):
                    sum_ = 0

                    # sum_(j=1) ^ (i - 1)(a_i_j * x_j(iter))
                    for j in range(0, i):
                        sum_ = sum_ + (self.A[i, j] * self.x[j])

                    # sum_(j=i + 1) ^ (n)(A_i_j * x_j(iter - 1))
                    for j in range(i+1, self.n):
                        sum_ = sum_ + (self.A[i, j] * self.x_old[j])

                    # x_i(iter) = (1-omega) * x_i(iter-1) + ((omega/a_i_i) * (b_i - s1 -s2)
                    self.x[i] = ((1-self.omega)*self.x_old[i]) + ((self.omega / self.A[i, i]) * (self.b[i] - sum_))

                # stopping criteria: (||x(iter) - x(iter-1)|| / ||x(iter)||) < TOL
                res = np.linalg.norm((self.x - self.x_old)) / np.linalg.norm(self.x_old)
                if res < self.TOL:
                    break

                # update the iteration counter and x(iter-1)
                self.x_old = self.x.copy()
                iter += 1

            # time taken for convergence
            time_taken = ti.time() - start_time

            # add the computed information to the Data class
            Data.add_data('x', self.x)
            Data.add_data('time_taken', time_taken)
            Data.add_data('iterations', iter-1)
            Data.add_data('residual', res)
