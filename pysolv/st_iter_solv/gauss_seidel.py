"""Gauss-Seidel Solver"""


# import the necessary packages
from data import *
import numpy as np

class GaussSeidelSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system"""

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self.init_x0()

        # call the gauss seidel solver
        self.solve()

    def init_x0(self):
        """Check if an initial value exists for the solution vector or initialize the solution vector"""

        try:
            self.x0
        except AttributeError:
            self.x0 = np.ones(self.n)

        # initialize the solution vectors
        self.x = np.empty(self.n)   # solution vector at the i'th iteration
        self.x_old = np.empty(self.n) # solution vector at the i-1'th iteration

    def solve(self):
        """Serial python implementation of Gauss Seidel solver."""

        # currently cannot handle non square systems in Gauss-seidel
        # ToDo: Implement generalized Gauss-seidel for non-square systems based on M.Saha's paper
        if self.m != self.n:
            msg = 'pysolv does not include a jacobi solver for non-square systems yet'
            raise ValueError(msg)

        # Gauss-Seidel iterations
        # x_i(iter) = (1/a_i_i)(b_i - sum_(j=1)^(i-1)(a_i_j * x_j(iter)) - sum_(j=i+1)^(n)(A_i_j * x_j(iter-1)))
        else:
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
                #u Update the iteration counter and x(iter-1)
                self.x_old = self.x.copy()
                iter += 1

            #add the solution to the Data class
            Data.add_data('x', self.x)

