"""Minimal residual algorithm.

Minimal residual algorithm for solving a set of linear equations is a type of 1D project scheme where the search space
and the left space spans only one vector i.e,

K -> span{v} ; L -> span{w}

Minimal residual algorithm can be used when the coefficient matrix (A) is not symmetric but positive definite.

For the MinRes algorithm, the search space spans {r} where as the left space spans {Ar}. The MinRes method minimizes the
following function

f(x) = norm2(b-A.x)**2
"""

# import the necessary packages
from pysolv.data import Data
import numpy as np


class MinResSolve(Data):

    def __init__(self):
        """Initialize the class with necessary inheritance and the solve the linear system.
        """

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self._init_x0()

        # call the conjugate gradient solver
        self._solve()

    def _init_x0(self):
        """Check if an initial value exists for the solution vector or initialize the solution vector
        """

        # initialize the solution vector
        try:
            self.x = self.x0
        except AttributeError:
            self.x = np.ones(self.n)

    def _py_solve(self):
        """Minimal residual iterations
        """

        # initialize the iteration counter
        iter = 0

        # initial residual and search direction
        # dir = res = b - A.x_0
        res = self.b - self.A.dot(self.x)
        dir_ = self.A.dot(res)

        # continue the iterations till itermax is reached
        while iter < self.itermax:
            # compute the step size, alpha
            alpha = (dir_.dot(res)) / (dir_.dot(dir_))

            # update the solution
            self.x = self.x + (alpha*res)

            # update the residual
            res = res - (alpha*dir_)

            # check for convergence
            if np.linalg.norm(res) <= self.tol:
                break

            # update the search direction
            dir_ = self.A.dot(res)

            # update the iteration counter
            iter += 1

        print(iter)
        print(np.linalg.norm(res))

    def _solve(self):
        """Serial implementation of the steepest descent iterative solver
        """

        if self.m != self.n:
            msg = 'Steepest descent solver cannot handle non-square systems'
            raise ValueError(msg)

        # steepest descent iterations
        else:
            # call the fortran sd solver
            #f_sdsolve.solve(self.A, self.b, self.x, self.itermax, self.tol)
            self._py_solve()

            # add the solution to the Data class
            Data.add_data('x', self.x)