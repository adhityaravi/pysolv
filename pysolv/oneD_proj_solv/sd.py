"""Steepest descent algorithm.

Steepest descent algorithm for solving a set of linear equations is a type of 1D project scheme where the search space
and the left space spans only one vector i.e,

K -> span{v} ; L -> span{w}

Steepest descent algorithm can be used when the coefficient matrix (A) is a SPD (symmetric positive definite) matrix.

At each iteration of the SD algorithm, the residual r is taken as v and w. And hence, the steepest descent algorithm
minimizes the following function

f(x) = (||x - x*||_A)**2 = (A.(x - x*), (x - x*))

Note: This solver could possibly have a worse convergence compared to even stationary iterative schemes like the SOR
"""

# import the necessary pacakges
from pysolv.data import Data
from pysolv.f_lib import f_sdsolve
import numpy as np


class SDSolve(Data):

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

    def _solve(self):
        """Serial implementation of the steepest descent iterative solver
        """

        if self.m != self.n:
            msg = 'Steepest descent solver cannot handle non-square systems'
            raise ValueError(msg)

        # steepest descent iterations
        else:
            # call the fortran sd solver
            f_sdsolve.solve(self.A, self.b, self.x, self.itermax, self.tol)

            # add the solution to the Data class
            Data.add_data('x', self.x)
