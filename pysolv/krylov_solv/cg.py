"""Conjugate Gradient linear solver.

"""

# import the necessary packages
from data import *
import numpy as np
import time as ti


class CGSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
           ToDO: Deprecate Conjugate Gradient _solve() and rather call SOR with omega=1
        """

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self._init_x0()

        # call the Conjugate Gradient solver
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
        self.x = self.x0

    def _solve(self):
        """Serial python implementation of Conjugate Gradient solver."""

        # currently cannot handle non square systems in Conjugate Gradient
        # ToDo: Implement generalized Conjugate Gradient for non-square systems based on M.Saha's paper
        if self.m != self.n:
            msg = 'pysolv does not include a Conjugate Gradient solver for non-square systems yet'
            raise ValueError(msg)

        # Conjugate Gradient iterations
        else:
            # time the solver
            start_time = ti.time()

            # Initialize the iteration counter and solution vector
            iter = 0

            # continue the iteration till the iteration counter reaches the maximum count if convergence is not obtained
            d = self.b - np.dot(self.A, self.x)
            r = d.copy()
            alpha_num = np.dot(r, r)
            print('Beginning CG')
            #while iter < 1:
            while iter < self.ITERMAX:

                ad = np.dot(self.A, d)
                alpha = alpha_num/np.dot(ad, d)
                self.x = self.x + alpha*d
                r = r - alpha*ad
                alpha_num_next = np.dot(r, r)
                beta = alpha_num_next/alpha_num
                d = r + beta*d
                alpha_num = alpha_num_next

                # stopping criteria: (||x(iter) - x(iter-1)|| / ||x(iter)||) < TOL
                # ToDo: Implement a pysolv native function to compute the residual
                res = np.linalg.norm((self.x - self.x_old)) / np.linalg.norm(self.x_old)

                self.x_old = self.x.copy()
                if res < self.TOL:
                    break
                # Update the iteration counter and x(iter-1)
                if iter % 10 == 0:
                    print('iter = ', iter)
                iter += 1

            # time taken for convergence
            time_taken = ti.time() - start_time
            print(np.dot(self.A, self.x))
            # add the solution to the Data class
            Data.add_data('x', self.x)
            Data.add_data('time_taken', time_taken)
            Data.add_data('iterations', iter - 1)
            Data.add_data('residual', res)