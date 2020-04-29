"""Symmetric successive over-relaxation linear solver

"""

# import the necessary packages
from data import *
import numpy as np
import time as ti


class SSORSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
        """

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self.init_x0()

        # check for relaxation parameter
        self.check_omega()

        # call the gauss seidel solver
        self.solve()

    def init_x0(self):
        """Check if an initial value exists for the solution vector or initialize the solution vector
        """

        try:
            self.x0
        except AttributeError:
            self.x0 = np.ones(self.n)

        # initialize the solution vectors
        self.x = np.empty(self.n)   # solution vector at the i'th iteration
        self.x_old = self.x0  # solution vector at the i-1'th iteration

    def check_omega(self):
        """Check if the user has prescribed a value for the relaxation parameter. If not, prescribe a value of 1.7 based
           on Dr. Rosics's lecture (Institute of Scientific Computing, TU Braunschweig):
           [https://www.tu-braunschweig.de/en/wire/teaching/previous-terms/winter-2016-17]
        """

        try:
            self.omega
        except AttributeError:
            self.omega = 1.7

    def solve(self):
        """Serial python implementation of successive over-relaxation solver based on Dr. Rosic's lectures (Institute of
           Scientific Computing, TU Braunschweig):
           [https://www.tu-braunschweig.de/en/wire/teaching/previous-terms/winter-2016-17]
        """