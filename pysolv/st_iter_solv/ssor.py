"""Symmetric successive over-relaxation (SSOR) linear solver

SSOR linear solver is a modification of the successive over-relaxation (SOR) scheme. SSOR scheme combines two SOR
sweeps (a forward sweep and a backward sweep). Check sor.py to understand the working of the SOR scheme.

The convergence of the SSOR scheme is generally worse than that of the SOR scheme. The iteration matrix (C) in SSOR is
similar to a symmetric matrix because of the nature of the scheme. Hence, the motivation behind the SSOR scheme is to
use it as a pre-conditioner for other iterative schemes with symmetric matrices.

For more information on SSOR scheme, check the following literature

[1]. https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
[2]. https://mathworld.wolfram.com/SymmetricSuccessiveOverrelaxationMethod.html
[3]. http://netlib.org/linalg/html_templates/node17.html
"""

# import the necessary packages
from pysolv.data import *
import numpy as np
import time as ti


class SSORSolve(Data):

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
        self.x_half = self.x0  # half solution from the first (forward) SOR sweep

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
        """Serial python implementation of successive over-relaxation solver based on Yousef Saad's book: Iterative
           method for sparse linear systems [https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf]
        """

        # currently the sor solver of pysolv cannot handle non-square systems
        if self.m != self.n:
            msg = 'pysolv does not include a ssor solver for non-square systems yet'
            raise ValueError(msg)

        # SSOR iteration
        # x(iter) = B1.B2.x(iter-1) + omega(2-omega)(1/(D + omega.U))(D)(1/(D + omega.L))b
        # B1 = (1/(D + omega.U))(-omega.L + (1-omega).D) --> Backward sweep
        # B2 = (1/(D + omega.L))(-omega.U + (1-omega).D) --> Forward sweep
        else:
            # time the solver
            start_time = ti.time()

            # initialize iteration counter
            iter = 0

            # continue the iteration till the iteration counter reaches the maximum count if convergence is not
            # obtained
            while iter < self.ITERMAX:
                # Forward sweep
                for i in range(0, self.n):
                    sum_ = 0

                    for j in range(0, i):
                        sum_ = sum_ + (self.A[i, j] * self.x_half[j])

                    for j in range(i+1, self.n):
                        sum_ = sum_ + (self.A[i, j] * self.x_old[j])

                    sum_ = (1/self.A[i, i])*(self.b[i] - sum_)
                    self.x_half[i] = (self.omega * sum_) + ((1-self.omega) * self.x_old[i])

                # Backward sweep
                for i in range(self.n-1, -1, -1):
                    sum_ = 0

                    for j in range(0, i):
                        sum_ = sum_ + (self.A[i, j] * self.x_half[j])

                    for j in range(i+1, self.n):
                        sum_ = sum_ + (self.A[i, j] * self.x[j])

                    sum_ = (1/self.A[i][i])*(self.b[i] - sum_)
                    self.x[i] = (self.omega * sum_) + ((1 - self.omega) * self.x_half[i])

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
            Data.add_data('iterations', iter - 1)
            Data.add_data('residual', res)
