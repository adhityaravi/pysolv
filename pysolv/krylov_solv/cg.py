"""Conjugate Gradient linear solver.
"""

# import the necessary packages
from pysolv import solve
from pysolv.data import Data
from pysolv.f_lib import f_cgsolve
import numpy as np


class CGSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
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
        """Conjugate gradient iterations
        """

        # initialize the iteration counter
        iter = 0

        # initial residual and search direction
        # dir = res = b - A.x_0
        res = self.b - self.A.dot(self.x)
        dir = res.copy()

        while iter < self.itermax:
            # compute the step size: alpha = res.res / dir.(A.dir)
            rdr = res.dot(res)
            Add = self.A.dot(dir)
            ddAdd = dir.dot(Add)
            alpha = rdr / ddAdd

            # compute the new solution
            self.x = self.x + (alpha * dir)

            # update the residual
            res_nxt = res - (alpha * Add)
            # check for convergence
            if np.linalg.norm(res_nxt) < self.tol:
                break

            # compute beta coefficient (from Gram-Schmidt orthogonalization)
            # beta = res_nxt.res_nxt / res.res
            beta = (res_nxt.dot(res_nxt)) / rdr

            # update the search direction
            dir = res_nxt + (beta * dir)
            res = res_nxt.copy()

            # update iteration counter
            iter += 1

    def _solve(self):
        """Serial implementation of Conjugate Gradient solver.
        """

        if self.m != self.n:
            msg = 'pysolv does not include a CG solver for non-square systems yet'
            raise ValueError(msg)

        # conjugate gradient iterations
        else:
            if self.n <= 400:
                # call the fortran cg solver
                f_cgsolve.init(self.A, self.b, self.tol, self.itermax)
                f_cgsolve.solve(self.x)
            else:
                # call the python cg solver
                self._py_solve()

            # add the solution to the Data class
            Data.add_data('x', self.x)
