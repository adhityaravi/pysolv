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

This solver also allows the user to use the SOR solver symmetrically --> Symmetric SOR.

SSOR linear solver is a modification of the successive over-relaxation (SOR) scheme. SSOR scheme combines two SOR
sweeps (a forward sweep and a backward sweep).

The convergence of the SSOR scheme is generally worse than that of the SOR scheme. The iteration matrix (C) in SSOR is
similar to a symmetric matrix because of the nature of the scheme. Hence, the motivation behind the SSOR scheme is to
use it as a pre-conditioner for other iterative schemes with symmetric matrices.

For more information on SSOR scheme, check the following literature

[1]. https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
[2]. https://mathworld.wolfram.com/SymmetricSuccessiveOverrelaxationMethod.html
[3]. http://netlib.org/linalg/html_templates/node17.html
"""

# import the necessary packages
from pysolv.f_lib import f_sorsolve
from pysolv.data import *
import numpy as np
import re


class SORSolve(Data):

    def __init__(self, symmetric_solve=False):
        """Initialize the class with the necessary inheritance and solve the linear system
        """

        # flag to switch between sor and symmetric sor
        self.symmetric_solve = symmetric_solve

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self._init_x0()

        # check for relaxation parameter and prepare omega related parameters for the fortran module
        self.adaptive_omega = 0
        self.h = 2
        self._map_omega()

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

    def _map_omega(self):
        """Check if the user has prescribed a value for the relaxation parameter. If not, compute omega adaptively based
           on Armijo conditions (Miyatake et al.)
        """

        # check if a relaxation parameter is already prescribed by the user
        try:
            # check if the omega should be computed adaptively
            if type(self.omega) == str:
                # raise error if adaptive omega technique asked for is not available
                if re.sub('\W+', '', self.omega).casefold() not in list(self.ADAPTIVE_OMEGA_FLAVOR.keys()):
                    msg = 'Please select from available adaptive omega techniques: {}'. \
                        format(list(self.ADAPTIVE_OMEGA_FLAVOR.keys()))
                    raise ValueError(msg)

                # convert the adaptive omega technique name (alphanumeric) to fortran recognized number
                else:
                    key = re.sub('\W+', '', self.omega).casefold()
                    self.adaptive_omega = self.ADAPTIVE_OMEGA_FLAVOR[key]
                    self.omega = 1

        # default relaxation parameter
        except AttributeError:
            self.adaptive_omega = 2
            self.omega = 1

        # default number of iterations after which relaxation parameter will be updated
        try:
            self.omega_update_frequency
        except AttributeError:
            self.omega_update_frequency = 20

    def _solve(self):
        """Serial implementation of successive over-relaxation solver based on Dr. Rosic's lectures (Institute of
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
            # call the sor fortran wrapper
            f_sorsolve.init(self.A, self.b, self.x0, self.itermax, self.tol, self.omega, self.h, self.adaptive_omega,
                            self.omega_update_frequency, self.c1, self.c2, self.lambda1, self.lambda2, self.rho1)

            # SSOR
            if self.symmetric_solve:
                f_sorsolve.s_solve(self.x)
            # SOR
            else:
                f_sorsolve.solve(self.x)

            # add the computed information to the Data class
            Data.add_data('x', self.x)
